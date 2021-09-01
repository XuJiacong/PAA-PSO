function [gbest,gbestval,hist_v,hist_p,FES]= DMS_PSO_func(fhd,Dimension,group_num,Particle_Number,iter_max,VRmin,VRmax,varargin)

    rand('state',sum(100*clock));
    group_ps=Particle_Number;
    ps=group_num*group_ps;
    D=Dimension;
    me = iter_max;
    L_FES=3*D;
    L_num=ceil(0.25.*group_num);
    cc=[1.49445 1.49445];
    iwt=0.729;

    if length(VRmin)==1
        VRmin=repmat(VRmin,1,D);
        VRmax=repmat(VRmax,1,D);
    end
    VRmin=repmat(VRmin,ps,1);
    VRmax=repmat(VRmax,ps,1);

    mv=0.2*(VRmax-VRmin);
    Vmin=-mv;
    Vmax=-Vmin;
    
    pos=VRmin+(VRmax-VRmin).*rand(ps,D);
    e=feval(fhd,pos',varargin{:});
    fitcount=ps;
    vel=Vmin+2.*Vmax.*rand(ps,D);
    pbest=pos;
    pbestval=e; 
    
    group_id = zeros(group_num,group_ps);
    pos_group = zeros(ps,1);
    gbestval = zeros(group_num,1);
    gbest = zeros(group_num, D);
    hist_v = zeros(me+1,1);
    hist_p = zeros(me+1,D);
    
    for i=1:group_num
        group_id(i,:)=((i-1)*group_ps+1):i*group_ps;
        pos_group(group_id(i,:))=i;
        %initialize the gbest and the gbest's fitness value for each group
        [gbestval(i),gbestid]=min(pbestval(group_id(i,:)));
        gbest(i,:)=pbest(group_id(i,gbestid),:);
    end

    threshold = ceil(0.95*me);
    for i = 1:threshold

        for k=1:ps
            aa=cc(1).*rand(1,D).*(pbest(k,:)-pos(k,:))+cc(2).*rand(1,D).*(gbest(pos_group(k),:)-pos(k,:));
            vel(k,:)=iwt.*vel(k,:)+aa; 
            vel(k,:)=(vel(k,:)>mv(k,:)).*mv(k,:)+(vel(k,:)<=mv(k,:)).*vel(k,:); 
            vel(k,:)=(vel(k,:)<(-mv(k,:))).*(-mv(k,:))+(vel(k,:)>=(-mv(k,:))).*vel(k,:);
            pos1=pos(k,:)+vel(k,:); 
            keep_d=rand(1,D)<0.5;
            pos(k,:)=keep_d.*pbest(k,:)+(1-keep_d).*pos1; 
            pos(k, :) = ((pos(k,:)<=VRmax(k,:))&(pos(k,:)>=VRmin(k,:))) .* pos(k,:) ...
                        + (pos(k,:)<VRmin(k,:)) .* (VRmin(k,:)+0.25.*(VRmax(k,:)-VRmin(k,:)).*rand(1,D)) ...
                        + (pos(k,:)>VRmax(k,:)) .* (VRmax(k,:)-0.25.*(VRmax(k,:)-VRmin(k,:)).*rand(1,D));    

            e(k)=feval(fhd,pos(k,:)',varargin{:});
            fitcount=fitcount+1;
            tmp=(pbestval(k)<e(k));
            temp=repmat(tmp,1,D);
            pbest(k,:)=temp.*pbest(k,:)+(1-temp).*pos(k,:);
            pbestval(k)=tmp.*pbestval(k)+(1-tmp).*e(k);%update the pbest
            if pbestval(k)<gbestval(pos_group(k))
                gbest(pos_group(k),:)=pbest(k,:);
                gbestval(pos_group(k))=pbestval(k);
            end

        end

        if mod(i,100)==0
            options = optimset('LargeScale','off','MaxFunEvals',L_FES,'Display','off');
            [~,tmpid]=sort(gbestval);
            for k=1:L_num
                [x,fval,~,output] = fminunc(fhd,gbest(tmpid(k),:)',options,varargin{:});
                fitcount=fitcount+output.funcCount;
                if fval<gbestval(tmpid(k))
                    [gbestval(tmpid(k)),gbestid]=min(pbestval(group_id(tmpid(k),:)));
                    pbest(group_id(tmpid(k),gbestid),:)=x';
                    pbestval(group_id(tmpid(k),gbestid))=fval;
                    gbest(tmpid(k),:)=x';
                    gbestval(tmpid(k))=fval;
                end
            end
        end

        if mod(i,5)==0
            rc=randperm(ps);
            group_id=zeros(group_num,group_ps);
            gbestval = zeros(group_num,1);
            gbest = zeros(group_num, D);
            for k=1:group_num
                group_id(k,:)=rc(((k-1)*group_ps+1):k*group_ps);
                pos_group(group_id(k,:))=k;
                [gbestval(k),gbestid]=min(pbestval(group_id(k,:)));
                gbest(k,:)=pbest(group_id(k,gbestid),:);
            end
        end

        [~,bbb]=sort(gbestval);
        hist_v(i+1) = gbestval(bbb(1));
        hist_p(i+1,:) = gbest(bbb(1));

    end
    [~,tmpid]=sort(gbestval);
    gbest=gbest(tmpid(1),:);
    gbestval=gbestval(tmpid(1));

    options = optimset('LargeScale','off','MaxFunEvals',0.05*me*ps,'Display','off');
    [x,fval,~,output] = fminunc(fhd,gbest',options,varargin{:});
    fitcount=fitcount+output.funcCount;
    if fval<gbestval
        gbest=x';
        gbestval=fval;
    end
    [~,bbb]=sort(gbestval);
    hist_v(i+1) = gbestval(bbb(1));
    hist_p(i+1,:) = gbest(bbb(1));
    
    for i = (threshold+1):me-1

        for k=1:ps
            aa=cc(1).*rand(1,D).*(pbest(k,:)-pos(k,:))+cc(2).*rand(1,D).*(gbest-pos(k,:));
            vel(k,:)=iwt.*vel(k,:)+aa; 
            vel(k,:)=(vel(k,:)>mv(k,:)).*mv(k,:)+(vel(k,:)<=mv(k,:)).*vel(k,:); 
            vel(k,:)=(vel(k,:)<(-mv(k,:))).*(-mv(k,:))+(vel(k,:)>=(-mv(k,:))).*vel(k,:);
            pos(k,:)=pos(k,:)+vel(k,:); 
            
            pos(k, :) = ((pos(k,:)<=VRmax(k,:))&(pos(k,:)>=VRmin(k,:))) .* pos(k,:) ...
                        + (pos(k,:)<VRmin(k,:)) .* (VRmin(k,:)+0.25.*(VRmax(k,:)-VRmin(k,:)).*rand(1,D)) ...
                        + (pos(k,:)>VRmax(k,:)) .* (VRmax(k,:)-0.25.*(VRmax(k,:)-VRmin(k,:)).*rand(1,D));    

            e(k)=feval(fhd,pos(k,:)',varargin{:});
            fitcount=fitcount+1;
            tmp=(pbestval(k)<e(k));
            temp=repmat(tmp,1,D);
            pbest(k,:)=temp.*pbest(k,:)+(1-temp).*pos(k,:);
            pbestval(k)=tmp.*pbestval(k)+(1-tmp).*e(k);%update the pbest
            if pbestval(k)<gbestval
            gbest=pbest(k,:);
            gbestval=pbestval(k);
            end
        end

        
        [~,bbb]=sort(gbestval);
        hist_v(i+1) = gbestval(bbb(1));
        hist_p(i+1,:) = gbest(bbb(1));
        
    end
    FES = fitcount;
end
