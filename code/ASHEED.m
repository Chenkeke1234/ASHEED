clear
FND=zeros(50,4);
HND=zeros(50,4);
Node.EnCCN=zeros(50,2);%Energy consumption of cluster heads and nodes
Node1.EnCCN=zeros(50,2);
Node2.EnCCN=zeros(50,2);
Node3.EnCCN=zeros(50,2);
Node.AvL=zeros(50,300);
Node1.AvL=zeros(50,300);
Node2.AvL=zeros(50,300);
Node3.AvL=zeros(50,300);
Node.Data=zeros(1,10);
Node1.Data=zeros(1,10);
Node2.Data=zeros(1,10);
Node3.Data=zeros(1,10);
for l=1:1
rmax=1000;
ETX=50*0.000000001;%传输能量，每bit
ERX=50*0.000000001;%接收能量，每bit
%Transmit Amplifier types
Efs=10*0.000000000001;%自由空间损耗，每bitclear
Emp=0.0013*0.000000000001;%多径衰落损耗，每bit
%Data Aggregation Energy
EDA=50*0.00000001;%融合能耗，每bit
EPC=50*0.00000001;%每轮的感知收集单位信息的能耗，每bit
%算出参数 do，节点通信半径
do=sqrt(Efs/Emp)/2;
RC0=100;
NodeNums=100; % the num of node 
n=NodeNums;
allive3=NodeNums;
d0=do;              %能耗模型参数
d1=do;
Elec=ERX; %      能耗参数
Eamp=Efs;  %
Eda=EDA;            %能量聚合能耗
Bx=0;  % The Postion of Baseation
By=0;
MaxInteral=20; % the leach simulate time                                            
Kbit=200; % 发送广播信息量
KB=4000;   %节点发送到簇头节点的单轮的数据量
Ers=EnRec(Elec,KB);                                                                %簇头节点对蔟内数据的数据融合系数
perdegree=pi/180*30;
InitEn=1.5;
NodeTranR=30;   %  the transit Radius
Pch=0.05;  % the desired percentage of cluster heads                          初始簇头比例  
P=Pch;
Tr=100;  % the time of round                                                   循环时间?
Gathingcoefficient=0.8;           %数据聚合系数
BandWitch = 1*10.^(6); %  Channel Bandwitch
Threshold=0;    % the threshold of node becoming a cluster-head                  
Cprob=0.04;    %初始簇头概率
NON_CH=0;%			// non cluster head
TENTATIVE_CH= 1; %			// tentative cluster head				
FINAL_CH= 2;%	// final cluster head
for i=1:NodeNums
Node.x(i)=(rand(1,1)*2-1)*n;  % the position of node               
Node.y(i)=(rand(1,1)*2-1)*n;
end
sym ClusterHeadNum ;          % 这个有意义吗？                   
TOS_LOCAL_ADDRESS = -1;       % TOS_LOCAL_ADDRESS  must <=0      判断自身为簇头情况      
Node.IsClusterHeads=linspace(0,0,NodeNums); % NON_CH,TENTATIVE_CH,FINAL_CH
Node.IsCovered=linspace(0,0,NodeNums);      % Have Been Covered by a cluster head 1:yes 0:No
Node.c=linspace(0,0,NodeNums);              % the Cluster head of node
Node.chcost=linspace(0,0,NodeNums);              % the Cluster head of node
Node.d=linspace(0,0,NodeNums);              % the distance between cluster head and node
Node.l=zeros(1,NodeNums)+Kbit;           % the length of node i transmit packet
Node.StateNode=ones(1,NodeNums);      % the State of all node 1: alive 0:dead
Node.Listothernode=zeros(NodeNums);         % if node is a cluster head,Listothernode save the id of node belong to this cluster       这个注意点对不对？
Node.csize=linspace(0,0,NodeNums);          % cluser size ,each cluster node num
Node.Nbr=zeros(NodeNums);                   % neighbor of node
Node.NumNbr=linspace(0,0,NodeNums);         % the neighbor's num of node
Node.DistNbr=linspace(0,0,NodeNums);         % the neighbor's dist of node                 
Node.CHprob=zeros(1,NodeNums); 
Node.InitCHprob=zeros(1,NodeNums);
Node.tent_CH=zeros(1,NodeNums)+NON_CH; 
Node.tent_CH_Cost=ones(1,NodeNums)+9999;
Node.IsaddDummyRound=linspace(0,0,NodeNums);                      
Node.n_finalCH=linspace(0,0,NodeNums);                                 
Node.ListfinalCH=zeros(NodeNums);                                      
Node.ListfinalCH_Cost=zeros(NodeNums)+9999;
Node.n_tentCH=linspace(0,0,NodeNums);                             
Node.ListtentCH=zeros(NodeNums);
Node.ListtentCH_Cost=zeros(NodeNums)+9999;
Node.my_finalCH=linspace(0,0,NodeNums);
Node.my_tentCH=linspace(0,0,NodeNums);
Node.my_final_CH_Cost=ones(1,NodeNums)+9999;
Node.Isstop=ones(1,NodeNums);    
AvEn=zeros(1,MaxInteral);
Avcost=zeros(1,MaxInteral);
AliveNode=zeros(1,MaxInteral);
MaxEn=zeros(1,NodeNums)+InitEn;  %最大能量值
Node.EnNode=zeros(1,NodeNums);
Node.EnNode=zeros(1,NodeNums)+InitEn;
Node1.x=Node.x;  % the position of node               
Node1.y=Node.y;
Node1.IsClusterHeads=linspace(0,0,NodeNums); % NON_CH,TENTATIVE_CH,FINAL_CH
Node1.IsCovered=linspace(0,0,NodeNums);      % Have Been Covered by a cluster head 1:yes 0:No
Node1.c=linspace(0,0,NodeNums);              % the Cluster head of node
Node1.chcost=linspace(0,0,NodeNums);              % the Cluster head of node
Node1.d=linspace(0,0,NodeNums);              % the distance between cluster head and node
Node1.l=zeros(1,NodeNums)+Kbit;           % the length of node i transmit packet
Node1.StateNode=ones(1,NodeNums);      % the State of all node 1: alive 0:dead
Node1.Listothernode=zeros(NodeNums);         % if node is a cluster head,Listothernode save the id of node belong to this cluster       这个注意点对不对？
Node1.csize=linspace(0,0,NodeNums);          % cluser size ,each cluster node num
Node1.Nbr=zeros(NodeNums);                   % neighbor of node
Node1.NumNbr=linspace(0,0,NodeNums);         % the neighbor's num of node
Node1.DistNbr=linspace(0,0,NodeNums);         % the neighbor's dist of node                 
Node1.CHprob=zeros(1,NodeNums); 
Node1.InitCHprob=zeros(1,NodeNums);
Node1.tent_CH=zeros(1,NodeNums)+NON_CH; 
Node1.tent_CH_Cost=ones(1,NodeNums)+9999;
Node1.IsaddDummyRound=linspace(0,0,NodeNums);                      
Node1.n_finalCH=linspace(0,0,NodeNums);                                 
Node1.ListfinalCH=zeros(NodeNums);                                      
Node1.ListfinalCH_Cost=zeros(NodeNums)+9999;
Node1.n_tentCH=linspace(0,0,NodeNums);                             
Node1.ListtentCH=zeros(NodeNums);
Node1.ListtentCH_Cost=zeros(NodeNums)+9999;
Node1.my_finalCH=linspace(0,0,NodeNums);
Node1.my_tentCH=linspace(0,0,NodeNums);
Node1.my_final_CH_Cost=ones(1,NodeNums)+9999;
Node1.Isstop=ones(1,NodeNums);    
MaxEn1=zeros(1,NodeNums)+InitEn;  %最大能量值
Node1.EnNode=zeros(1,NodeNums);
Node1.EnNode=zeros(1,NodeNums)+InitEn;
Node2.x=Node.x;  % the position of node               
Node2.y=Node.y;
Node2.IsClusterHeads=linspace(0,0,NodeNums); % NON_CH,TENTATIVE_CH,FINAL_CH
Node2.IsCovered=linspace(0,0,NodeNums);      % Have Been Covered by a cluster head 1:yes 0:No
Node2.c=linspace(0,0,NodeNums);              % the Cluster head of node
Node2.chcost=linspace(0,0,NodeNums);              % the Cluster head of node
Node2.d=linspace(0,0,NodeNums);              % the distance between cluster head and node
Node2.l=zeros(1,NodeNums)+Kbit;           % the length of node i transmit packet
Node2.StateNode=ones(1,NodeNums);      % the State of all node 1: alive 0:dead
Node2.Listothernode=zeros(NodeNums);         % if node is a cluster head,Listothernode save the id of node belong to this cluster       这个注意点对不对？
Node2.csize=linspace(0,0,NodeNums);          % cluser size ,each cluster node num
Node2.Nbr=zeros(NodeNums);                   % neighbor of node
Node2.NumNbr=linspace(0,0,NodeNums);         % the neighbor's num of node
Node2.DistNbr=linspace(0,0,NodeNums);         % the neighbor's dist of node                 
Node2.CHprob=zeros(1,NodeNums); 
Node2.InitCHprob=zeros(1,NodeNums);
Node2.tent_CH=zeros(1,NodeNums)+NON_CH; 
Node2.tent_CH_Cost=ones(1,NodeNums)+9999;
Node2.IsaddDummyRound=linspace(0,0,NodeNums);                      
Node2.n_finalCH=linspace(0,0,NodeNums);                                 
Node2.ListfinalCH=zeros(NodeNums);                                      
Node2.ListfinalCH_Cost=zeros(NodeNums)+9999;
Node2.n_tentCH=linspace(0,0,NodeNums);                             
Node2.ListtentCH=zeros(NodeNums);
Node2.ListtentCH_Cost=zeros(NodeNums)+9999;
Node2.my_finalCH=linspace(0,0,NodeNums);
Node2.my_tentCH=linspace(0,0,NodeNums);
Node2.my_final_CH_Cost=ones(1,NodeNums)+9999;
Node2.Isstop=ones(1,NodeNums);    
MaxEn2=zeros(1,NodeNums)+InitEn;  %最大能量值
Node2.EnNode=zeros(1,NodeNums);
Node2.EnNode=zeros(1,NodeNums)+InitEn;
Node3.x=Node.x;
Node3.y=Node.y;
Node3.EnNode=zeros(1,NodeNums);
Node3.EnNode=zeros(1,NodeNums)+InitEn;
Node3.G=zeros(1,n);
Node3.D=zeros(1,n);
Node3.IsClusterHeads=linspace(0,0,NodeNums);
Node3.NumNbr=linspace(0,0,NodeNums);
Node3.Cid=zeros(1,n);
% FND=zeros(50,4);
% HND=zeros(50,4);
Flag_FND=zeros(1,4);
Flag_HND=zeros(1,4);
data=0;
data1=0;
data2=0;
data3=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%建立邻居节点表
Node3.NB=zeros(n,n);
for i=1:1:n
    for j=1:1:n
        d=sqrt((Node3.x(i)-(Node3.x(j)))^2 + (Node3.y(i)-(Node3.y(j)))^2);
        if(d<=NodeTranR && d~=0 && Node3.EnNode(j)>0)%节点间距离少于感知半径皆为邻居节点
            Node3.NB(i,j)=j;
        end
    end
    dist=Node3.x(i)*Node3.x(i)+Node3.y(i)*Node3.y(i);
    EntranPCH=EnTran(Elec,Eamp,Kbit,dist);
    Node3.EnNode(i)=Node3.EnNode(i)-EntranPCH;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LEACH-N算法
for r=1:1:rmax     %该 for 循环将下面的所有程序包括在内，直到最后一 end 才结束循环
    r
  %每过一个轮转周期(本程序为20次)使各节点的S（i）.G参数（该参数用于后面的簇选举，
  % 在该轮转周期内已当选过簇头的节点不能再当选）恢复为零，即50个周期，每周期20轮
  if(mod(r,round(1/P))==0)%round四舍五入取整函数, mod(x,y)为求余函数，y为除数，符号与除数相同
    for i=1:1:n
        Node3.G(i)=0;
    end
  end 
  Et3=0;
  for i=1:1:n
      if(Node3.EnNode(i)>0)
        Et3=Node3.EnNode(i)+Et3;%El3为剩余总能量
      end
  end
%(2)死亡节点检查模块
dead3=0;
for i=1:1:n
    %检查有无死亡节点
    if (Node3.EnNode(i)<=0)
        dead3=dead3+1;
        if(dead3==1)
            if(Flag_FND(4)==0)
                FND(l,4)=r;
                Flag_FND(4)=1;
            end
        end
        if(dead3==0.5*NodeNums)
            if(Flag_HND(4)==0)
                HND(l,4)=r;
                Flag_HND(4)=1;
            end
        end
    end
end
STATISTICS.DEAD3(r)=dead3;%死亡节点总数
STATISTICS.ALLIVE3(r)=allive3-dead3;%存活节点数
m3=allive3-dead3;
Ea3=Et3/(allive3-dead3);%Et为总能量，Et*(1-r/rmax)为每一轮剩余的总能量，Ea为剩余总平均能量
fnd3=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%求邻节点个数
 for i=1:NodeNums
     if(Node3.EnNode(i)>0)
         count =0 ;
        for j=1:NodeNums  
            if(j~=i && Node3.EnNode(j)>0) 
            dist = ((Node3.x(i)-Node3.x(j)).^2)+((Node3.y(i)-Node3.y(j)).^2);  % the distance.^2
                   if (dist<NodeTranR^2)           % 考虑cost设置不同的时候需要将dist的值写入进去
                       count=count+1;
                       Node3.Nbr(i,count)=j;       
                   end    
             end
             if j== NodeNums 
                    Node3.NumNbr(i) = count ;          %邻居节点总个数
             end  
        end 
     end
 end 
%(4)簇头选举模块
countCHs3=0;
cluster3=1;
C3=[];
for i=1:1:n
 if Ea3>0
 p(i)=P*Node3.NumNbr(i)/n;
 if(Node3.EnNode(i)>0)
   temp_rand=rand;     
   if ((Node3.G(i))<=0)  
       %簇头的选举，当选的簇头会把各种相关信存入下面程序所给定的变量中
        if(temp_rand<= (p(i)/(1-p(i)*mod(r,round(1/p(i))))*(Node3.EnNode(i)/InitEn)))%p(i)<1,大概率
            %p(i)/(1-p(i)*mod(r,round(1/p(i)))))为阈值计算公式
            countCHs3=countCHs3+1;
            Node3.IsClusterHeads(i)=FINAL_CH;
            Node3.Cid(i)=0;
            Node3.G(i)=round(1/p(i))-1;%p(i)越小，S3(i).G越大，大于0时不能竞选簇头
            C3(cluster3)=i;
            cluster3=cluster3+1;
            dist=Node3.x(i)*Node3.x(i)+Node3.y(i)*Node3.y(i);
            EntranPCH=EnTran(Elec,Eamp,Kbit,dist);
            Node3.EnNode(i)=Node3.EnNode(i)-EntranPCH;
            Node3.EnCCN(l,1)=Node3.EnCCN(l,1)+EntranPCH;
        end      
    end
    % S3(i).G=S3(i).G-1;     
  end 
 end
end
STATISTICS.COUNTCHS5(r)=countCHs3;%簇头个数

%(5)簇内成员选择簇头模块(即簇的形成模块)
%簇内成员对簇头的选择（即簇的形成）算法
LB5=zeros(1,n);
for i=1:1:n
   if ( Node3.IsClusterHeads(i)==FINAL_CH && Node3.EnNode(i)>0 )
       min_dis=Inf;%iInf表示正无穷大
       for c=1:1:cluster3-1
           temp=min(min_dis,sqrt((Node3.x(i)-Node3.x(C3(c)))^2 + (Node3.y(i)-Node3.y(C3(c)))^2));
           if ( temp<min_dis )
               min_dis=temp;
               Node3.Cid(i)=C3(c);%保存普通所连接的簇头的编号
               LB5(Node3.Cid(i))=LB5(Node3.Cid(i))+1;%计算簇内节点个数
           end
       end
%        if(min_dis<S5(i).D)
%            LB5(S5(i).Cid)=LB5(S5(i).Cid)+1;%计算簇内节点个数
%        else
%            S5(i).Cid='sink';
%        end
   end
end
%%%%%%计算建立连接的能耗
for i=1:n
    if(Node3.IsClusterHeads(i)==FINAL_CH)
       EnRecP=EnRec(Elec,Kbit*LB5(i));
       Node3.EnNode(i)=Node3.EnNode(i)-EnRecP;
       Node3.EnCCN(l,1)=Node3.EnCCN(l,1)+EnRecP;
    else
        if(Node3.Cid(i)~=0)
            dist=Node3.x(i)*Node3.x(Node3.Cid(i))+Node3.y(i)*Node3.y(Node3.Cid(i));
            EntranPCH=EnTran(Elec,Eamp,Kbit,dist);
            Node3.EnNode(i)=Node3.EnNode(i)-EntranPCH;
        else
            dist=Node3.x(i)*Node3.x(i)+Node3.y(i)*Node3.y(i);
            EntranPCH=EnTran(Elec,Eamp,Kbit,dist);
            Node3.EnNode(i)=Node3.EnNode(i)-EntranPCH;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%能耗模型
 for i=1:1:n
     if(Node3.EnNode(i)>0)
         if(Node3.Cid(i)~=0 )
             d=sqrt((Node3.x(i)-Node3.x(Node3.Cid(i)))^2 + (Node3.y(i)-Node3.y(Node3.Cid(i)))^2);
            if (d>do)
                Node3.EnNode(i)=Node3.EnNode(i)-((EPC+ETX)*(4000)+(ERX+EDA)*(4000)*LB5(i) + Emp*4000*( d * d * d * d));
            end
            if (d<=do)
                Node3.EnNode(i)=Node3.EnNode(i)-( (EPC+ETX)*(4000)+(ERX+EDA)*(4000)*LB5(i) + Efs*4000*( d * d)); 
            end
         else
              d=sqrt(Node3.x(i)^2 + Node3.y(i)^2);
              if (d>do)
                    Node3.EnNode(i)=Node3.EnNode(i)-( ETX*(4000)+(ERX+EDA)*(4000)*LB5(i) + Emp*4000*( d * d * d * d)); 
              end
                if (d<=do)
                    Node3.EnNode(i)=Node3.EnNode(i)-( ETX*(4000)+(ERX+EDA)*(4000)*LB5(i) + Efs*4000*( d * d)); 
                end
         end   
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%计算簇头能耗%%%%%%%%%%%%             
         if(Node3.IsClusterHeads(i)==FINAL_CH)
             if(Node3.Cid(i)~=0)
                d=sqrt((Node3.x(i)-Node3.x(Node3.Cid(i)))^2 + (Node3.y(i)-Node3.y(Node3.Cid(i))^2));
                if(d>do)
                   Node3.EnCCN(l,1)=Node3.EnCCN(l,1)+((EPC+ETX)*(4000)+(ERX+EDA)*(4000)*LB5(i)+Emp*4000*(d * d * d * d));
                end
                if (d<=do)
                    Node3.EnCCN(l,1)=Node3.EnCCN(l,1)+((EPC+ETX)*(4000)+(ERX+EDA)*(4000)*LB5(i)+Efs*4000*(d * d)); 
                end
             else
                d=sqrt((Node3.x(i))^2 + (Node3.y(i))^2);
                if(d>do)
                   Node3.EnCCN(l,1)=Node3.EnCCN(l,1)+((EPC+ETX)*(4000)+(ERX+EDA)*(4000)*LB5(i)+Emp*4000*(d * d * d * d));
                end
                if (d<=do)
                    Node3.EnCCN(l,1)=Node3.EnCCN(l,1)+((EPC+ETX)*(4000)+(ERX+EDA)*(4000)*LB5(i)+Efs*4000*(d * d)); 
                end
             end
         end
     end 
 end
Et3=0;
for i=1:1:n
  if(Node3.EnNode(i)>0)
    Et3=Node3.EnNode(i)+Et3;%El3为剩余总能量
  end
end
Node3.EnCCN(l,2)=Et3;
for i=1:1:n
  if(Node3.EnNode(i)<=0 && Node3.AvL(l,i)==0)
    Node3.AvL(l,i)=r;
    if(dead3>=0.5*NodeNums)
        for j=1:NodeNums
            if(Node3.AvL(l,j)==0)
                Node3.AvL(l,j)=r;
            end
        end
    end
  end
end
for i=1:n
    if(Node3.EnNode(i)>0 && dead3<=0.5*n)
        data3=data3+1;
    end
end
end
Node3.Data(1,l)=data3;


%%%%%%%%%%%%%%%%%%%%%%%%%ASHEED
%%%计算全集群竞争半径
DQRC1=[];
Node1.QRC=zeros(1,NodeNums); 
for i=1:NodeNums
    %%%%%计算簇头与基站的最大值与最小值
    for j=1:NodeNums
        DQRC1(j)=sqrt(Node1.x(j)*Node1.x(j)+Node1.y(j)*Node1.y(j));
    end
end
%%%%%%计算竞争半径
for j=1:NodeNums
    Node1.QRC(j)=(1-(max(DQRC1)-sqrt(Node1.x(j)*Node1.x(j)+Node1.y(j)*Node1.y(j)))/(3*(max(DQRC1)-min(DQRC1))))*RC0;
    dist=Node1.x(j)*Node1.x(j)+Node1.y(j)*Node1.y(j);
    EntranPCH=EnTran(Elec,Eamp,Kbit,dist);
    Node1.EnNode(j)=Node1.EnNode(j)-EntranPCH; %将竞争半径发送到基站，进一步确定FINAL_CH
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%求邻节点个数
 for i=1:NodeNums
     count =0 ;
    for j=1:NodeNums  
        if(j~=i) 
        dist = ((Node1.x(i)-Node1.x(j)).^2)+((Node1.y(i)-Node1.y(j)).^2);  % the distance.^2
               if (dist < Node1.QRC(i)^2)           % 考虑cost设置不同的时候需要将dist的值写入进去
                   count=count+1;
                   Node1.Nbr(i,count)=j;       
               end    
         end
         if j== NodeNums 
                Node1.NumNbr(i) = count ;          %邻居节点总个数
         end  
    end 
 end                                               %初始化时求出节点的邻居节点及邻节点个数 
 sym iteration;
for i=1:NodeNums
    Node1.CHprob(i)=Cprob*((Node1.EnNode(i))/MaxEn1(i));
end
ClusterHeadNum=0;
iteration=0;
dead1=0;
allive=NodeNums;
%选簇
while sum(Node1.Isstop)~=0
    iteration=iteration+1;
    for i =1:NodeNums             %此时进行选簇迭代过程，由于仿真每次都有节点id为1到NodeNums进行，导致id号靠前的节点其当选簇头概率较大，能耗较多。改进时可将节点id每次随机化
        if Node1.Isstop(i)==1                      
          if Node1.CHprob(i)<1         
             if Node1.tent_CH(i)==NON_CH
              if rand(1,1)<Node1.CHprob(i)
                 Node1.IsClusterHeads(i)=TENTATIVE_CH; 
                 Node1.tent_CH(i)=TOS_LOCAL_ADDRESS;
                 Node1.tent_CH_Cost(i)=Node1.NumNbr(i);                                  %cost值设置为节点的邻节点数目？即以度作为其归蔟度量值
              end
              %elseif  Node.tent_CH(i)==TOS_LOCAL_ADDRESS
             end
           Node1.CHprob(i)=Node1.CHprob(i).*2;
          else            %即CHprob等于1时的情况
             for j=1:NodeNums   %Node.n_finalCH(i)
                if Node1.ListfinalCH(i,j) ~=0
                 if Node1.my_final_CH_Cost(i) > Node1.ListfinalCH_Cost(i,j) 
                     Node1.my_finalCH(i)= Node1.ListfinalCH(i,j);                           %进行归蔟行动
                     Node1.my_final_CH_Cost(i)=Node1.ListfinalCH_Cost(i,j);
                 end
                end  
             end
              % choose cluster head
              Node1.Isstop(i)=0;   % diedai end   until CHprob==1
              if Node1.my_finalCH(i) ~= NON_CH
                 Node1.IsClusterHeads(i)= NON_CH;                                                                      %这个为什么？
                 Node1.c(i)=Node1.my_finalCH(i);
                 Node1.chcost=Node1.my_final_CH_Cost(i);
                 %join the cluster 
                 dist =((Node1.x(i)-Node1.x(curentnbr)).^2)+((Node1.y(i)-Node1.y(curentnbr)).^2); % ((Node.x(i)-Node.x(curentnbr)).^2)+((Node.y(i)-Node.y(curentnbr)).^2);  % the distance.^2                这个值大概用来广播信息时用的
                 EntranPCH=EnTran(Elec,Eamp,Kbit,dist) ;                                                          %大概是传播能量损耗值计算 函数
                 Node1.d(i)=((Node1.x(i)-Node1.x(Node1.c(i))).^2)+((Node1.y(i)-Node1.y(Node1.c(i))).^2);  % the distance.^2
                 Node1.EnNode(i)=Node1.EnNode(i)-EntranPCH;                 %预先判断其能量是否足以进行这次传输
                 if Node1.EnNode(i) <= 0
                    Node1.StateNode(i)=0;
                    Node1.Isstop(i)=0;
                    Node1.EnNode(i)=0;
                 end
                EnRecP=EnRec(Elec,Kbit);                                          %这个是簇头接收数据后的能耗计算函数
                Node1.EnNode(Node1.c(i))=Node1.EnNode(Node1.c(i))-EnRecP;                
                if Node1.EnNode(Node1.c(i)) <= 0
                    Node1.StateNode(Node1.c(i))=0;
                    Node1.Isstop(Node1.c(i))=0;
                    Node1.EnNode(Node1.c(i))=0; 
                else                   
                    Node1.csize(Node1.c(i))=Node1.csize(Node1.c(i))+1;  % cluster size add one
                end
             else
                 Node1.IsClusterHeads(i)= FINAL_CH;
                 Node1.my_finalCH(i)=TOS_LOCAL_ADDRESS;
                 Node1.c(i)=TOS_LOCAL_ADDRESS;
                 Node1.my_final_CH_Cost(i)= Node1.NumNbr(i);%computeDegree(i);
                 Node1.chcost=Node1.my_final_CH_Cost(i);
                 Node1.d(i)=((Node1.x(i)-Bx).^2)+((Node1.y(i)-By).^2);  % the distance.^2  %簇头节点到sink节点的距离
                 ClusterHeadNum=ClusterHeadNum+1;
              end
          end    
       end
      % compute consume energy
    if Node1.IsClusterHeads(i) == TENTATIVE_CH  % & Node.tent_CH(curentnbr)==TOS_LOCAL_ADDRESS
        dist =NodeTranR.^2; % ((Node.x(i)-Node.x(curentnbr)).^2)+((Node.y(i)-Node.y(curentnbr)).^2);  % the distance.^2                    大概需要改正
           EntranPCH=EnTran(Elec,Eamp,Kbit,dist) ;                                 %即减去广播信息能耗
           Node1.EnNode(i)=Node1.EnNode(i)-EntranPCH;
           if Node1.EnNode(i) <= 0
                Node1.StateNode(i)=0;
                Node1.Isstop(i)=0;
                Node1.EnNode(i) =0;
           end
        for j=1:Node1.NumNbr(i)
            curentnbr = Node1.Nbr(i,j);                %处在临时簇头状态的蔟发送广播信息 求其邻居节点
             EnRecP=EnRec(Elec,Kbit);
           Node1.EnNode(curentnbr)=Node1.EnNode(curentnbr)-EnRecP;                  %这点是否需要加上在临时簇头节点能量值大于0的情况下
           if Node1.EnNode(curentnbr) > 0
              if (Node1.ListtentCH(curentnbr,i)==0)
                   Node1.n_tentCH(curentnbr)=Node1.n_tentCH(curentnbr)+1;    
              end 
              Node1.ListtentCH(curentnbr,i) = i;
              Node1.ListtentCH_Cost(curentnbr,i)=Node1.NumNbr(i);%Node.computeDegree(i);

              % if Node.tent_CH(curentnbr)~=TOS_LOCAL_ADDRESS &(Node.tent_CH(curentnbr)==NON_CH | Node.tent_CH_Cost(i)< Node.tent_CH_Cost(curentnbr)  | ((Node.tent_CH_Cost(i)== Node.tent_CH_Cost(curentnbr) ) & i < Node.tent_CH(curentnbr)))
                if(Node1.tent_CH(curentnbr)==NON_CH || Node1.tent_CH_Cost(i)< Node1.tent_CH_Cost(curentnbr)  || ((Node1.tent_CH_Cost(i)== Node1.tent_CH_Cost(curentnbr) ) && i < Node1.tent_CH(curentnbr)))
                    Node1.tent_CH_Cost(curentnbr)=Node1.tent_CH_Cost(i);
                    Node1.tent_CH(curentnbr)=i;
                end     % if(Node.tent_CH(curentnbr)==NON_CH | Node.tent_CH_Cost(i)< Node.tent_CH_Cost(curentnbr)
           else
                 Node1.StateNode(curentnbr)=0;
                 Node1.Isstop(curentnbr)=0;
                 Node1.EnNode(curentnbr)=0; 
           end  
           end
            elseif  Node1.IsClusterHeads(i) == FINAL_CH
                dist = NodeTranR.^2; %((Node.x(i)-Node.x(curentnbr)).^2)+((Node.y(i)-Node.y(curentnbr)).^2);  % the distance.^2
                   EntranPCH=EnTran(Elec,Eamp,Kbit,dist);
                   Node1.EnNode(i)=Node1.EnNode(i)-EntranPCH;
                   if Node1.EnNode(i) <= 0
                        Node1.StateNode(i)=0;
                        Node1.Isstop(i)=0;
                        Node1.EnNode(i)=0;
                   end  
            for j=1:Node1.NumNbr(i)
                curentnbr = Node1.Nbr(i,j);
                EnRecP=EnRec(Elec,Kbit);
                Node1.EnNode(curentnbr)=Node1.EnNode(curentnbr)-EnRecP;
               if Node1.EnNode(curentnbr) > 0
                   if (Node1.ListfinalCH(curentnbr,i)==0)                                   %这个是否跟上面一样错误
                       Node1.n_finalCH(curentnbr)=Node1.n_finalCH(curentnbr)+1;    
                   end
                  Node1.ListfinalCH(curentnbr,i)=i;
                  Node1.ListfinalCH_Cost(curentnbr,i)=Node1.NumNbr(i);%Node.computeDegree(i);
               else
                 Node1.EnNode(curentnbr)=0; 
                 Node1.StateNode(curentnbr)=0;
                 Node1.Isstop(curentnbr)=0;
               end    
            end
        end
      end  
 end  
%%%%%%%%%聚类
DRC1=[];
Rc1=[];
Rc1=find(Node1.IsClusterHeads==FINAL_CH);
RCSET1=zeros(NodeNums,NodeNums);%保存簇头簇内成员节点ID
for i=1:NodeNums
    %%%%%计算簇头与基站的最大值与最小值
    for j=1:NodeNums
        if Node1.IsClusterHeads(i) == FINAL_CH
            DRC1(j)=sqrt(Node1.x(j)*Node1.x(j)+Node1.y(j)*Node1.y(j));
        end
    end
end
%%%%计算簇头竞争半径
for j=1:NodeNums
    if (Node1.IsClusterHeads(j) == FINAL_CH)
        Node1.RC(j)=(1-(max(DRC1)-sqrt(Node1.x(j)*Node1.x(j)+Node1.y(j)*Node1.y(j)))/(3*(max(DRC1)-min(DRC1))))*RC0;
        dist=Node1.x(j)*Node1.x(j)+Node1.y(j)*Node1.y(j);
        EntranPCH=EnTran(Elec,Eamp,Kbit,dist);
        Node1.EnNode(j)=Node1.EnNode(j)-EntranPCH; %将竞争半径发送到基站，进一步确定FINAL_CH
        Node1.EnCCN(l,1)=Node1.EnCCN(l,1)+EntranPCH;
    end
end
%%%%%利用竞争半径减少簇头节点数量
for j=1:NodeNums
    if Node1.IsClusterHeads(j) == FINAL_CH
        for z=1:1:NodeNums
            if(Node1.IsClusterHeads(z) == FINAL_CH && j~=z)
                df=sqrt((Node1.x(j)-Node1.x(z))*(Node1.x(j)-Node1.x(z)) + (Node1.y(j)-Node1.y(z))*(Node1.y(j)-Node1.y(z)));
                if (df>Node1.RC(j) && df>Node1.RC(z))
                    Node1.IsClusterHeads(j)=FINAL_CH;
                    Node1.IsClusterHeads(z)=FINAL_CH;
                else
                    Node1.IsClusterHeads(z)=NON_CH;
                    Node1.IsClusterHeads(j)=NON_CH;
                end
            end
        end
    EnRecP=EnRec(Elec,Kbit);
    Node1.EnNode(j)=Node1.EnNode(j)-EnRecP;
    Node1.EnCCN(l,1)=Node1.EnCCN(l,1)+EnRecP;
    end
end

%%%%%%%%%%统计簇头竞争半径内的节点
for j=1:NodeNums
    if (Node1.IsClusterHeads(j) == FINAL_CH)
        for z=1:NodeNums
            Djz=sqrt((Node1.x(j)-Node1.x(z))*(Node1.x(j)-Node1.x(z)) + (Node1.y(j)-Node1.y(z))*(Node1.y(j)-Node1.y(z)));
            if(Node1.RC(j)>=Djz && j~=z)
               RCSET1(j,z)=z;  %%保存簇头RC内的节点ID
            end
        end
    end
end

%%%%%聚类
CH1=[];
CSET1=zeros(NodeNums,NodeNums);%记录簇集成员节点ID
CHcount1=length(find(Node1.IsClusterHeads==2));
CH1=find(Node1.IsClusterHeads==2);
LBRU1=zeros(1,NodeNums);
%%%%%%%%建立连接
for i=1:NodeNums
    if (sqrt(Node1.x(i)^2 + Node1.y(i)^2)<=do)
        Node1.Cid(i)=0;%与基站连接
    else
        if(Node1.IsClusterHeads(i)==NON_CH && Node1.EnNode(i)>0)
           min_dis=Inf;%iInf表示正无穷大
           dist=Node1.x(i)*Node1.x(i)+Node1.y(i)*Node1.y(i);
           EntranPCH=EnTran(Elec,Eamp,Kbit,dist);
           Node1.EnNode(i)=Node1.EnNode(i)-EntranPCH;
           for c=1:1:CHcount1
               temp=min(min_dis,sqrt((Node1.x(CH1(c))-Node1.x(i))^2 + (Node1.y(CH1(c))-Node1.y(i))^2));
               if (temp<min_dis)
                   min_dis=temp;
                   Node1.Cid(i)=CH1(c);%保存普通所连接的簇头的编号
               end
           end
           CSET1(Node1.Cid(i),i)=i;%保存对应簇头簇集内成员节点
           LBRU1(Node1.Cid(i))=LBRU1(Node1.Cid(i))+1;
        end
    end
end
%%%%%%%%%%%簇头在RC中选择最优传输路径
for i=1:NodeNums
    if(Node1.IsClusterHeads(i)==FINAL_CH && Node1.EnNode(i)>0)
       min_dis=Inf;%iInf表示正无穷大
       for c=1:1:length(Rc1)
           temp=min(min_dis,(sqrt((Node1.x(Rc1(c))-Node1.x(i))^2 + (Node1.y(Rc1(c))-Node1.y(i))^2)+Node1.x(Rc1(c))^2+Node1.y(Rc1(c))^2));
           if (temp<min_dis)
               min_dis=temp;
               Node1.Cid(i)=Rc1(c);%保存普通所连接的簇头的编号
           end
       end
       LBRU1(Node1.Cid(i))=LBRU1(Node1.Cid(i))+1;
    end
end

%%%%%%%%计算普通节点加入簇集时的簇头接收能耗
for i=1:NodeNums
    Kbit=200;
    if(Node1.IsClusterHeads(i) == FINAL_CH)
        Kbit=Kbit*LBRU1(i);
        EnRecP=EnRec(Elec,Kbit);
        Node1.EnNode(i)=Node1.EnNode(i)-EnRecP;
        Node1.EnCCN(l,1)=Node1.EnCCN(l,1)+EnRecP;
    end
end

for r=1:1:rmax
dead1=0;
r
Es=0;
for i=1:NodeNums
    if(Node1.EnNode(i)>0)
        Es=Node1.EnNode(i)+Es;
    else
        dead1=dead1+1;
    end
end
dead1=0;
deadCH1=0;
for i=1:1:NodeNums
    if (Node1.EnNode(i)<=0)
        dead1=dead1+1;
        if(dead1==1)
            if(Flag_FND(2)==0)
                FND(l,2)=r;
                Flag_FND(2)=1;
            end
        end
        if(dead1==0.5*NodeNums)
            if(Flag_HND(2)==0)
                HND(l,2)=r;
                Flag_HND(2)=1;
            end
        end
        if(Node1.IsClusterHeads(i)==FINAL_CH)
           deadCH1=1;
        end
    end
end
STATISTICS.DEAD1(r)=dead1;%死亡节点总数
STATISTICS.ALLIVE1(r)=allive-dead1;%存活节点数

DQRC1=[];
Node1.QRC=zeros(1,NodeNums); 
for i=1:NodeNums
    %%%%%计算簇头与基站的最大值与最小值
    for j=1:NodeNums
        DQRC1(j)=sqrt(Node1.x(j)*Node1.x(j)+Node1.y(j)*Node1.y(j));
    end
end
%%%%%%计算竞争半径
for j=1:NodeNums
    Node1.QRC(j)=(1-(max(DQRC1)-sqrt(Node1.x(j)*Node1.x(j)+Node1.y(j)*Node1.y(j)))/(3*(max(DQRC1)-min(DQRC1))))*RC0;
    dist=Node1.x(j)*Node1.x(j)+Node1.y(j)*Node1.y(j);
    EntranPCH=EnTran(Elec,Eamp,Kbit,dist);
    Node1.EnNode(j)=Node1.EnNode(j)-EntranPCH; %将竞争半径发送到基站，进一步确定FINAL_CH
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%更新邻节点个数
 for i=1:NodeNums
     if(Node1.EnNode(j)>0)
         count =0;
        for j=1:NodeNums  
            if(j~=i && Node1.EnNode(j)>0) 
            dist = ((Node1.x(i)-Node1.x(j)).^2)+((Node1.y(i)-Node1.y(j)).^2);  % the distance.^2
                   if (dist < Node1.QRC(i)^2)            % 考虑cost设置不同的时候需要将dist的值写入进去
                       count=count+1;
                       Node1.Nbr(i,count)=j;       
                   end    
             end
             if j== NodeNums 
                    Node1.NumNbr(i) = count ;          %邻居节点总个数
             end  
        end
     end
 end 

m1=0;
if(dead1>0)
    m1=STATISTICS.DEAD1(r)-STATISTICS.DEAD1(r-1);
end
HalfCH1=zeros(1,NodeNums);%%%保存能耗少于网络平均能量的簇头ID
HalfECH1=0;
for i=1:1:NodeNums
    if(Node1.IsClusterHeads(i)==FINAL_CH && Node1.EnNode(i)<(Es/(allive-dead1)))
       if(Node1.EnNode(i)<=0)
           HalfECH1=0;
       else
           HalfCH1(i)=i; %保存簇头剩余能量大于0而少于平均能耗的ID
       end
    end
end

HC=[];
HC=find(HalfCH1>0);%簇头节点ID
if(m1>0)
while sum(Node1.Isstop)~=0
    iteration=iteration+1;
    for i =1:NodeNums             %此时进行选簇迭代过程，由于仿真每次都有节点id为1到NodeNums进行，导致id号靠前的节点其当选簇头概率较大，能耗较多。改进时可将节点id每次随机化
        if Node1.Isstop(i)==1                      
          if Node1.CHprob(i)<1         
             if Node1.tent_CH(i)==NON_CH
              if rand(1,1)<Node1.CHprob(i)
                 Node1.IsClusterHeads(i)=TENTATIVE_CH; 
                 Node1.tent_CH(i)=TOS_LOCAL_ADDRESS;
                 Node1.tent_CH_Cost(i)=Node1.NumNbr(i);                                  %cost值设置为节点的邻节点数目？即以度作为其归蔟度量值
              end
              %elseif  Node.tent_CH(i)==TOS_LOCAL_ADDRESS
             end
           Node1.CHprob(i)=Node1.CHprob(i).*2;
          else            %即CHprob等于1时的情况
             for j=1:NodeNums   %Node.n_finalCH(i)
                if Node1.ListfinalCH(i,j) ~=0
                 if Node1.my_final_CH_Cost(i) > Node1.ListfinalCH_Cost(i,j) 
                     Node1.my_finalCH(i)= Node1.ListfinalCH(i,j);                           %进行归蔟行动
                     Node1.my_final_CH_Cost(i)=Node1.ListfinalCH_Cost(i,j);
                 end
                end  
             end
              % choose cluster head
            Node1.Isstop(i)=0;   % diedai end   until CHprob==1
            if Node1.my_finalCH(i) ~= NON_CH
                Node1.IsClusterHeads(i)= NON_CH;                                                                      %这个为什么？
                Node1.c(i)=Node1.my_finalCH(i);
                Node1.chcost=Node1.my_final_CH_Cost(i);
                %join the cluster 
                dist =((Node1.x(i)-Node1.x(curentnbr)).^2)+((Node1.y(i)-Node1.y(curentnbr)).^2); % ((Node.x(i)-Node.x(curentnbr)).^2)+((Node.y(i)-Node.y(curentnbr)).^2);  % the distance.^2                这个值大概用来广播信息时用的
                EntranPCH=EnTran(Elec,Eamp,Kbit,dist) ;                                                          %大概是传播能量损耗值计算 函数
                Node1.d(i)=((Node1.x(i)-Node1.x(Node1.c(i))).^2)+((Node1.y(i)-Node1.y(Node1.c(i))).^2);  % the distance.^2
                Node1.EnNode(i)=Node1.EnNode(i)-EntranPCH;                 %预先判断其能量是否足以进行这次传输
                if Node1.EnNode(i) <= 0
                    Node1.StateNode(i)=0;
                    Node1.Isstop(i)=0;
                    Node1.EnNode(i)=0;
                end
                EnRecP=EnRec(Elec,Kbit);                                          %这个是簇头接收数据后的能耗计算函数
                Node1.EnNode(Node1.c(i))=Node1.EnNode(Node1.c(i))-EnRecP;                
                if Node1.EnNode(Node1.c(i)) <= 0
                    Node1.StateNode(Node1.c(i))=0;
                    Node1.Isstop(Node1.c(i))=0;
                    Node1.EnNode(Node1.c(i))=0; 
                else                   
                    Node1.csize(Node1.c(i))=Node1.csize(Node1.c(i))+1;  % cluster size add one
                end
            else
                Node1.IsClusterHeads(i)= FINAL_CH;
                Node1.my_finalCH(i)=TOS_LOCAL_ADDRESS;
                Node1.c(i)=TOS_LOCAL_ADDRESS;
                Node1.my_final_CH_Cost(i)= Node1.NumNbr(i);%computeDegree(i);
                Node1.chcost=Node1.my_final_CH_Cost(i);
                Node1.d(i)=((Node1.x(i)-Bx).^2)+((Node1.y(i)-By).^2);  % the distance.^2  %簇头节点到sink节点的距离
                ClusterHeadNum=ClusterHeadNum+1;
            end
          end    
       end
      % compute consume energy
    if Node1.IsClusterHeads(i) == TENTATIVE_CH  % & Node.tent_CH(curentnbr)==TOS_LOCAL_ADDRESS
        dist =NodeTranR.^2; % ((Node.x(i)-Node.x(curentnbr)).^2)+((Node.y(i)-Node.y(curentnbr)).^2);  % the distance.^2                    大概需要改正
           EntranPCH=EnTran(Elec,Eamp,Kbit,dist) ;                                 %即减去广播信息能耗
           Node1.EnNode(i)=Node1.EnNode(i)-EntranPCH;
           if Node1.EnNode(i) <= 0
                Node1.StateNode(i)=0;
                Node1.Isstop(i)=0;
                Node1.EnNode(i) =0;
           end
        for j=1:Node1.NumNbr(i)
            curentnbr = Node1.Nbr(i,j);                %处在临时簇头状态的蔟发送广播信息 求其邻居节点
             EnRecP=EnRec(Elec,Kbit);
           Node1.EnNode(curentnbr)=Node1.EnNode(curentnbr)-EnRecP;                  %这点是否需要加上在临时簇头节点能量值大于0的情况下
           if Node1.EnNode(curentnbr) > 0
              if (Node1.ListtentCH(curentnbr,i)==0)
                   Node1.n_tentCH(curentnbr)=Node1.n_tentCH(curentnbr)+1;    
              end 
              Node1.ListtentCH(curentnbr,i) = i;
              Node1.ListtentCH_Cost(curentnbr,i)=Node1.NumNbr(i);%Node.computeDegree(i);

              % if Node.tent_CH(curentnbr)~=TOS_LOCAL_ADDRESS &(Node.tent_CH(curentnbr)==NON_CH | Node.tent_CH_Cost(i)< Node.tent_CH_Cost(curentnbr)  | ((Node.tent_CH_Cost(i)== Node.tent_CH_Cost(curentnbr) ) & i < Node.tent_CH(curentnbr)))
                if(Node1.tent_CH(curentnbr)==NON_CH || Node1.tent_CH_Cost(i)< Node1.tent_CH_Cost(curentnbr)  || ((Node1.tent_CH_Cost(i)== Node1.tent_CH_Cost(curentnbr) ) && i < Node1.tent_CH(curentnbr)))
                    Node1.tent_CH_Cost(curentnbr)=Node1.tent_CH_Cost(i);
                    Node1.tent_CH(curentnbr)=i;
                end     % if(Node.tent_CH(curentnbr)==NON_CH | Node.tent_CH_Cost(i)< Node.tent_CH_Cost(curentnbr)
           else
                 Node1.StateNode(curentnbr)=0;
                 Node1.Isstop(curentnbr)=0;
                 Node1.EnNode(curentnbr)=0; 
           end  
       end
            elseif  Node1.IsClusterHeads(i) == FINAL_CH
                dist = NodeTranR.^2; %((Node.x(i)-Node.x(curentnbr)).^2)+((Node.y(i)-Node.y(curentnbr)).^2);  % the distance.^2
                   EntranPCH=EnTran(Elec,Eamp,Kbit,dist);
                   Node1.EnNode(i)=Node1.EnNode(i)-EntranPCH;
                   if Node1.EnNode(i) <= 0
                        Node1.StateNode(i)=0;
                        Node1.Isstop(i)=0;
                        Node1.EnNode(i)=0;
                   end  
            for j=1:Node1.NumNbr(i)
                curentnbr = Node1.Nbr(i,j);
                EnRecP=EnRec(Elec,Kbit);
                Node1.EnNode(curentnbr)=Node1.EnNode(curentnbr)-EnRecP;
               if Node1.EnNode(curentnbr) > 0
                   if (Node1.ListfinalCH(curentnbr,i)==0)                                   %这个是否跟上面一样错误
                       Node1.n_finalCH(curentnbr)=Node1.n_finalCH(curentnbr)+1;    
                   end
                  Node1.ListfinalCH(curentnbr,i)=i;
                  Node1.ListfinalCH_Cost(curentnbr,i)=Node1.NumNbr(i);%Node.computeDegree(i);
               else
                 Node1.EnNode(curentnbr)=0; 
                 Node1.StateNode(curentnbr)=0;
                 Node1.Isstop(curentnbr)=0;
               end    
            end
    end
     end  
end  
    %%%%%%%%%聚类
DRC1=[];
Rc1=[];
Rc1=find(Node1.IsClusterHeads==FINAL_CH);
for i=1:NodeNums
    %%%%%计算簇头与基站的最大值与最小值
    for j=1:NodeNums
        if Node1.IsClusterHeads(i) == FINAL_CH
            DRC1(j)=sqrt(Node1.x(j)*Node1.x(j)+Node1.y(j)*Node1.y(j));
        end
    end
end
    %%%%%%计算竞争半径
for j=1:NodeNums
    if (Node1.IsClusterHeads(i) == FINAL_CH)
        Node1.RC(j)=(1-(max(DRC1)-sqrt(Node1.x(j)*Node1.x(j)+Node1.y(j)*Node1.y(j)))/(3*(max(DRC1)-min(DRC1))))*RC0;
        dist=Node1.x(j)*Node1.x(j)+Node1.y(j)*Node1.y(j);
        EntranPCH=EnTran(Elec,Eamp,Kbit,dist);
        Node1.EnNode(j)=Node1.EnNode(j)-EntranPCH;
        Node1.EnCCN(l,1)=Node1.EnCCN(l,1)+EntranPCH;
    end
end
%%%%%利用竞争半径减少簇头节点数量
for j=1:NodeNums
    if Node1.IsClusterHeads(j) == FINAL_CH
        for z=1:1:NodeNums
            if(Node1.IsClusterHeads(z) == FINAL_CH && j~=z)
                df=sqrt((Node1.x(j)-Node1.x(z))*(Node1.x(j)-Node1.x(z)) + (Node1.y(j)-Node1.y(z))*(Node1.y(j)-Node1.y(z)));
                if (df>Node1.RC(j) && df>Node1.RC(z))
                    Node1.IsClusterHeads(j)=FINAL_CH;
                    Node1.IsClusterHeads(z)=FINAL_CH;
                else
                    Node1.IsClusterHeads(z)=NON_CH;
                    Node1.IsClusterHeads(j)=NON_CH;
                end
            end
        end
    EnRecP=EnRec(Elec,Kbit);
    Node1.EnNode(j)=Node1.EnNode(j)-EnRecP;
    Node1.EnCCN(l,1)=Node1.EnCCN(l,1)+EnRecP;
    end
end
CH1=[];
CSET1=zeros(NodeNums,NodeNums);%记录簇集成员节点ID
CHcount1=length(find(Node1.IsClusterHeads==2));
CH1=find(Node1.IsClusterHeads==2);
LBRU1=zeros(1,NodeNums);
%%%%%%%%建立连接
for i=1:NodeNums
    if (sqrt(Node1.x(i)^2 + Node1.y(i)^2)<=do)
        Node1.Cid(i)=0;%与基站连接
    else
        if(Node1.IsClusterHeads(i)==NON_CH && Node1.EnNode(i)>0)
           min_dis=Inf;%iInf表示正无穷大
           dist=Node1.x(i)*Node1.x(i)+Node1.y(i)*Node1.y(i);
           EntranPCH=EnTran(Elec,Eamp,Kbit,dist);
           Node1.EnNode(i)=Node1.EnNode(i)-EntranPCH;
           for c=1:1:CHcount1
               temp=min(min_dis,sqrt((Node1.x(CH1(c))-Node1.x(i))^2 + (Node1.y(CH1(c))-Node1.y(i))^2));
               if (temp<min_dis)
                   min_dis=temp;
                   Node1.Cid(i)=CH1(c);%保存普通所连接的簇头的编号
               end
           end   
           CSET1(Node1.Cid(i),i)=i;%保存对应簇头簇集内成员节点
           LBRU1(Node1.Cid(i))=LBRU1(Node1.Cid(i))+1;%计算簇内节点个数
        end
    end
end
%%%%%%%%%%%簇头在RC中选择最优传输路径
for i=1:NodeNums
    if(Node1.IsClusterHeads(i)==FINAL_CH && Node1.EnNode(i)>0)
       min_dis=Inf;%iInf表示正无穷大
       for c=1:1:length(Rc1)
           temp=min(min_dis,(sqrt((Node1.x(Rc1(c))-Node1.x(i))^2 + (Node1.y(Rc1(c))-Node1.y(i))^2)+Node1.x(Rc1(c))^2+Node1.y(Rc1(c))^2));
           if (temp<min_dis)
               min_dis=temp;
               Node1.Cid(i)=Rc1(c);%保存普通所连接的簇头的编号
           end
       end
       LBRU1(Node1.Cid(i))=LBRU1(Node1.Cid(i))+1;
    end
end
end
if(HalfECH1>0) 
for j=1:length(HC)
    s=[];
    t=1;
    ES=0;
    temp=Inf;
    SC=0;
    for z=1:NodeNums 
        if(CSET1(HC(j),z)~=0)
            s(t)=CSET1(HC(j),z);
            t=t+1;
        end
    end
    for z=1:length(s)%总能量
        ES=Node1.EnNode(s(z))+ES;
    end
    for z=1:length(s)
        dist=(Node1.x(s(z))-Node1.x(HC(j)))^2+(Node1.y(s(z))-Node1.y(HC(j)))^2;
        EntranPCH=EnTran(Elec,Eamp,Kbit,dist);
        for i=1:length(s)
            if(i~=z)
                dist=(Node1.x(s(z))-Node1.x(s(i)))^2+(Node1.y(s(i))-Node1.y(s(i)))^2;
                EntranPCH=EnTran(Elec,Eamp,Kbit,dist)+EntranPCH;
            end
        end
        if(EntranPCH<temp)
            temp=EntranPCH;
            SC=s(z);
        end
    end
    if(Node1.EnNode(SC)>ES/length(s))
        for z=1:length(s)
            if(s(z)~=SC)
                dist1=(Node1.x(s(z))-Node1.x(SC))^2+(Node1.y(s(z))-Node1.y(SC))^2;
                dist2=(Node1.x(s(z))-Node1.x(HC(j)))^2+(Node1.y(s(z))-Node1.y(HC(j)))^2;
                if(dist1<dist2)
                    Node1.Cid(s(z))=SC;
                    EntranPCH=EnTran(Elec,Eamp,Kbit,dist1);
                    Node1.EnNode(s(z))=Node1.EnNode(s(z))-EntranPCH;
                    LBRU1(HC(j))=LBRU1(HC(j))-1;
                    LBRU1(SC)=LBRU1(SC)+1;
                end
            end
        end
        EnRecP=EnRec(Elec,Kbit);
        Node1.EnNode(SC)=Node1.EnNode(SC)-EnRecP;
        Node1.EnCCN(l,1)=Node1.EnCCN(l,1)+EnRecP;
    end
end
end

for i=1:NodeNums
    Kbit=200;
    if(Node1.IsClusterHeads(i) == FINAL_CH)
        Kbit=Kbit*LBRU1(i);
        EnRecP=EnRec(Elec,Kbit);
        Node1.EnNode(i)=Node1.EnNode(i)-EnRecP;
        Node1.EnCCN(l,1)=Node1.EnCCN(l,1)+EnRecP;
    end
end
%%%%%%%%%%%%%%%%%%%%%能耗模型%%%%%%
for i=1:NodeNums
     if(Node1.EnNode(i)>0)
         if(Node1.IsClusterHeads(i)==FINAL_CH && Node1.EnNode(i)<(Es/(allive-dead1)))
            d=sqrt((Node1.x(i)-Node1.x(Node1.Cid(i)))^2 + (Node1.y(i)-Node1.y(Node1.Cid(i))^2));
            if(d>do)
               Node1.EnNode(i)=Node1.EnNode(i)-((EPC+ETX)*(4000)+Emp*4000*(d * d * d * d));%+(ERX+EDA)*(4000)*LBRU1(i);
            end
            if (d<=do)
                Node1.EnNode(i)=Node1.EnNode(i)-((EPC+ETX)*(4000)+Efs*4000*(d * d));%+(ERX+EDA)*(4000)*LBRU1(i); 
            end
         else
             if(Node1.Cid(i)~=0)
                 d=sqrt((Node1.x(i)-Node1.x(Node1.Cid(i)))^2 + (Node1.y(i)-Node1.y(Node1.Cid(i))^2));
                if(d>do)
                   Node1.EnNode(i)=Node1.EnNode(i)-((EPC+ETX)*(4000)+Emp*4000*(d * d * d * d))+(ERX+EDA)*(4000)*LBRU1(i);
                end
                if (d<=do)
                    Node1.EnNode(i)=Node1.EnNode(i)-((EPC+ETX)*(4000)+Efs*4000*(d * d))+(ERX+EDA)*(4000)*LBRU1(i);
                end
             else
                d=sqrt((Node1.x(i))^2 + (Node1.y(i))^2);
                if(d>do)
                   Node1.EnNode(i)=Node1.EnNode(i)-((EPC+ETX)*(4000)+(ERX+EDA)*(4000)*LBRU1(i)+Emp*4000*(d * d * d * d));
                end
                if (d<=do)
                    Node1.EnNode(i)=Node1.EnNode(i)-((EPC+ETX)*(4000)+(ERX+EDA)*(4000)*LBRU1(i)+Efs*4000*(d * d)); 
                end
             end
         end
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%计算簇头能耗%%%%%%%%%%%%             
         if(Node1.IsClusterHeads(i)==FINAL_CH)
             if(Node1.Cid(i)~=0)
                 if(Node1.EnNode(i)<(Es/(allive-dead1)))
                     d=sqrt((Node1.x(i)-Node1.x(Node1.Cid(i)))^2 + (Node1.y(i)-Node1.y(Node1.Cid(i))^2));
                    if(d>do)
                       Node1.EnCCN(l,1)=Node1.EnCCN(l,1)+((EPC+ETX)*(4000)+Emp*4000*(d * d * d * d));
                    end
                    if (d<=do)
                        Node1.EnCCN(l,1)=Node1.EnCCN(l,1)+((EPC+ETX)*(4000)+Efs*4000*(d * d)); 
                    end
                 else
                    d=sqrt((Node1.x(i)-Node1.x(Node1.Cid(i)))^2 + (Node1.y(i)-Node1.y(Node1.Cid(i))^2));
                    if(d>do)
                       Node1.EnCCN(l,1)=Node1.EnCCN(l,1)+((EPC+ETX)*(4000)+(ERX+EDA)*(4000)*LBRU1(i)+Emp*4000*(d * d * d * d));
                    end
                    if (d<=do)
                        Node1.EnCCN(l,1)=Node1.EnCCN(l,1)+((EPC+ETX)*(4000)+(ERX+EDA)*(4000)*LBRU1(i)+Efs*4000*(d * d)); 
                    end
                 end
             else
                d=sqrt((Node1.x(i))^2 + (Node1.y(i))^2);
                if(Node1.EnNode(i)<(Es/(allive-dead1)))
                    if(d>do)
                       Node1.EnCCN(l,1)=Node1.EnCCN(l,1)+((EPC+ETX)*(4000)+Emp*4000*(d * d * d * d));
                    end
                    if (d<=do)
                        Node1.EnCCN(l,1)=Node1.EnCCN(l,1)+((EPC+ETX)*(4000)+Efs*4000*(d * d)); 
                    end
                else
                    if(d>do)
                       Node1.EnCCN(l,1)=Node1.EnCCN(l,1)+((EPC+ETX)*(4000)+(ERX+EDA)*(4000)*LBRU1(i)+Emp*4000*(d * d * d * d));
                    end
                    if (d<=do)
                        Node1.EnCCN(l,1)=Node1.EnCCN(l,1)+((EPC+ETX)*(4000)+(ERX+EDA)*(4000)*LBRU1(i)+Efs*4000*(d * d)); 
                    end
                end
             end
         end
     end
end
Es=0;
for i=1:NodeNums
    if(Node1.EnNode(i)>0)
        Es=Node1.EnNode(i)+Es;
    end
end
Node1.EnCCN(l,2)=Es;
for i=1:1:n
  if(Node1.EnNode(i)<=0 && Node1.AvL(l,i)==0)
    Node1.AvL(l,i)=r;
    if(dead1>=0.5*NodeNums)
        for j=1:NodeNums
            if(Node1.AvL(l,j)==0)
                Node1.AvL(l,j)=r;
            end
        end
    end
  end
end
for i=1:n
    if(Node1.EnNode(i)>0 && dead1<=0.5*n)
        data1=data1+1;
    end
end
end
Node1.Data(1,l)=data1;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RUHEED
%求邻节点个数
 for i=1:NodeNums
     count =0 ;
    for j=1:NodeNums  
        if(j~=i) 
        dist = ((Node.x(i)-Node.x(j)).^2)+((Node.y(i)-Node.y(j)).^2);  % the distance.^2
               if dist < NodeTranR^2           %       考虑cost设置不同的时候需要将dist的值写入进去
                   count=count+1;
                   Node.Nbr(i,count)=j;       
               end    
         end
         if j== NodeNums 
                Node.NumNbr(i) = count ;          %邻居节点总个数
         end  
    end 
 end                                               %初始化时求出节点的邻居节点及邻节点个数
 
 sym iteration;
for i=1:NodeNums
    Node.CHprob(i)=Cprob*((Node.EnNode(i))/MaxEn(i));
end
ClusterHeadNum=0;
iteration=0;
sum1=sum(Node.EnNode);
Node.RC=zeros(1,NodeNums); 
dead=0;
allive=NodeNums;
%选簇
while sum(Node.Isstop)~=0
    iteration=iteration+1;
    for i =1:NodeNums             %此时进行选簇迭代过程，由于仿真每次都有节点id为1到NodeNums进行，导致id号靠前的节点其当选簇头概率较大，能耗较多。改进时可将节点id每次随机化
        if Node.Isstop(i)==1                      
          if Node.CHprob(i)<1         
             if Node.tent_CH(i)==NON_CH
              if rand(1,1)<Node.CHprob(i)
                 Node.IsClusterHeads(i)=TENTATIVE_CH; 
                 Node.tent_CH(i)=TOS_LOCAL_ADDRESS;
                 Node.tent_CH_Cost(i)=Node.NumNbr(i);                                  %cost值设置为节点的邻节点数目？即以度作为其归蔟度量值
              end
              %elseif  Node.tent_CH(i)==TOS_LOCAL_ADDRESS
             end
           Node.CHprob(i)=Node.CHprob(i).*2;
          else            %即CHprob等于1时的情况
             for j=1:NodeNums   %Node.n_finalCH(i)
                if Node.ListfinalCH(i,j) ~=0
                 if Node.my_final_CH_Cost(i) > Node.ListfinalCH_Cost(i,j) 
                     Node.my_finalCH(i)= Node.ListfinalCH(i,j);                           %进行归蔟行动
                     Node.my_final_CH_Cost(i)=Node.ListfinalCH_Cost(i,j);
                 end
                end  
             end
              % choose cluster head
              Node.Isstop(i)=0;   % diedai end   until CHprob==1
              if Node.my_finalCH(i) ~= NON_CH
                 Node.IsClusterHeads(i)= NON_CH;                                                                      %这个为什么？
                 Node.c(i)=Node.my_finalCH(i);
                 Node.chcost=Node.my_final_CH_Cost(i);
                 %join the cluster 
                 dist =((Node.x(i)-Node.x(curentnbr)).^2)+((Node.y(i)-Node.y(curentnbr)).^2); % ((Node.x(i)-Node.x(curentnbr)).^2)+((Node.y(i)-Node.y(curentnbr)).^2);  % the distance.^2                这个值大概用来广播信息时用的
                 EntranPCH=EnTran(Elec,Eamp,Kbit,dist) ;                                                          %大概是传播能量损耗值计算 函数
                 Node.d(i)=((Node.x(i)-Node.x(Node.c(i))).^2)+((Node.y(i)-Node.y(Node.c(i))).^2);  % the distance.^2
                 Node.EnNode(i)=Node.EnNode(i)-EntranPCH;                 %预先判断其能量是否足以进行这次传输
                 if Node.EnNode(i) <= 0
                    Node.StateNode(i)=0;
                    Node.Isstop(i)=0;
                    Node.EnNode(i)=0;
                 end
                EnRecP=EnRec(Elec,Kbit);                                          %这个是簇头接收数据后的能耗计算函数
                Node.EnNode(Node.c(i))=Node.EnNode(Node.c(i))-EnRecP;                
                if Node.EnNode(Node.c(i)) <= 0
                    Node.StateNode(Node.c(i))=0;
                    Node.Isstop(Node.c(i))=0;
                    Node.EnNode(Node.c(i))=0; 
                else                   
                    Node.csize(Node.c(i))=Node.csize(Node.c(i))+1;  % cluster size add one
                end
             else
                 Node.IsClusterHeads(i)= FINAL_CH;
                 Node.my_finalCH(i)=TOS_LOCAL_ADDRESS;
                 Node.c(i)=TOS_LOCAL_ADDRESS;
                 Node.my_final_CH_Cost(i)= Node.NumNbr(i);%computeDegree(i);
                 Node.chcost=Node.my_final_CH_Cost(i);
                 Node.d(i)=((Node.x(i)-Bx).^2)+((Node.y(i)-By).^2);  % the distance.^2  %簇头节点到sink节点的距离
                 ClusterHeadNum=ClusterHeadNum+1;
              end
          end    
       end
      % compute consume energy
    if Node.IsClusterHeads(i) == TENTATIVE_CH  % & Node.tent_CH(curentnbr)==TOS_LOCAL_ADDRESS
        dist =NodeTranR.^2; % ((Node.x(i)-Node.x(curentnbr)).^2)+((Node.y(i)-Node.y(curentnbr)).^2);  % the distance.^2                    大概需要改正
           EntranPCH=EnTran(Elec,Eamp,Kbit,dist) ;                                 %即减去广播信息能耗
           Node.EnNode(i)=Node.EnNode(i)-EntranPCH;
           if Node.EnNode(i) <= 0
                Node.StateNode(i)=0;
                Node.Isstop(i)=0;
                Node.EnNode(i) =0;
           end
        for j=1:Node.NumNbr(i)
            curentnbr = Node.Nbr(i,j);                %处在临时簇头状态的蔟发送广播信息 求其邻居节点
             EnRecP=EnRec(Elec,Kbit);
           Node.EnNode(curentnbr)=Node.EnNode(curentnbr)-EnRecP;                  %这点是否需要加上在临时簇头节点能量值大于0的情况下
           if Node.EnNode(curentnbr) > 0
              if (Node.ListtentCH(curentnbr,i)==0)
                   Node.n_tentCH(curentnbr)=Node.n_tentCH(curentnbr)+1;    
              end 
              Node.ListtentCH(curentnbr,i) = i;
              Node.ListtentCH_Cost(curentnbr,i)=Node.NumNbr(i);%Node.computeDegree(i);

              % if Node.tent_CH(curentnbr)~=TOS_LOCAL_ADDRESS &(Node.tent_CH(curentnbr)==NON_CH | Node.tent_CH_Cost(i)< Node.tent_CH_Cost(curentnbr)  | ((Node.tent_CH_Cost(i)== Node.tent_CH_Cost(curentnbr) ) & i < Node.tent_CH(curentnbr)))
                if(Node.tent_CH(curentnbr)==NON_CH || Node.tent_CH_Cost(i)< Node.tent_CH_Cost(curentnbr)  || ((Node.tent_CH_Cost(i)== Node.tent_CH_Cost(curentnbr) ) && i < Node.tent_CH(curentnbr)))
                    Node.tent_CH_Cost(curentnbr)=Node.tent_CH_Cost(i);
                    Node.tent_CH(curentnbr)=i;
                end     % if(Node.tent_CH(curentnbr)==NON_CH | Node.tent_CH_Cost(i)< Node.tent_CH_Cost(curentnbr)
           else
                 Node.StateNode(curentnbr)=0;
                 Node.Isstop(curentnbr)=0;
                 Node.EnNode(curentnbr)=0; 
           end  
           end
            elseif  Node.IsClusterHeads(i) == FINAL_CH
                dist = NodeTranR.^2; %((Node.x(i)-Node.x(curentnbr)).^2)+((Node.y(i)-Node.y(curentnbr)).^2);  % the distance.^2
                   EntranPCH=EnTran(Elec,Eamp,Kbit,dist);
                   Node.EnNode(i)=Node.EnNode(i)-EntranPCH;
                   if Node.EnNode(i) <= 0
                        Node.StateNode(i)=0;
                        Node.Isstop(i)=0;
                        Node.EnNode(i)=0;
                   end  
            for j=1:Node.NumNbr(i)
                curentnbr = Node.Nbr(i,j);
                EnRecP=EnRec(Elec,Kbit);
                Node.EnNode(curentnbr)=Node.EnNode(curentnbr)-EnRecP;
               if Node.EnNode(curentnbr) > 0
                   if (Node.ListfinalCH(curentnbr,i)==0)                                   %这个是否跟上面一样错误
                       Node.n_finalCH(curentnbr)=Node.n_finalCH(curentnbr)+1;    
                   end
                  Node.ListfinalCH(curentnbr,i)=i;
                  Node.ListfinalCH_Cost(curentnbr,i)=Node.NumNbr(i);%Node.computeDegree(i);
               else
                 Node.EnNode(curentnbr)=0; 
                 Node.StateNode(curentnbr)=0;
                 Node.Isstop(curentnbr)=0;
               end    
            end
        end
      end  
 end  
%%%%%%%%%聚类
DRC=[];
Rc=[];
Rc=find(Node.IsClusterHeads==FINAL_CH);
for i=1:NodeNums
    %%%%%计算簇头与基站距离的最大值与最小值
    for j=1:NodeNums
        if Node.IsClusterHeads(i) == FINAL_CH
            DRC(j)=sqrt(Node.x(j)*Node.x(j)+Node.y(j)*Node.y(j));
        end
    end
end
%%%%%%计算竞争半径
for j=1:NodeNums
    if (Node.IsClusterHeads(j) == FINAL_CH)
        Node.RC(j)=(1-(max(DRC)-sqrt(Node.x(j)*Node.x(j)+Node.y(j)*Node.y(j)))/(3*(max(DRC)-min(DRC))))*RC0;
        dist=Node.x(j)*Node.x(j)+Node.y(j)*Node.y(j);
        EntranPCH=EnTran(Elec,Eamp,Kbit,dist);
        Node.EnNode(j)=Node.EnNode(j)-EntranPCH; %将竞争半径发送到基站，进一步确定FINAL_CH
        Node.EnCCN(l,1)=Node.EnCCN(l,1)+EntranPCH;
    end
end

%%%%%%利用竞争半径减少簇头节点数量
for j=1:NodeNums
    if Node.IsClusterHeads(j) == FINAL_CH
        for z=1:1:NodeNums
            if(Node.IsClusterHeads(z) == FINAL_CH && j~=z)
                df=sqrt((Node.x(j)-Node.x(z))*(Node.x(j)-Node.x(z)) + (Node.y(j)-Node.y(z))*(Node.y(j)-Node.y(z)));
                if (df>Node.RC(j) && df>Node.RC(z))
                    Node.IsClusterHeads(j)=FINAL_CH;
                    Node.IsClusterHeads(z)=FINAL_CH;
                else
                    Node.IsClusterHeads(z)=NON_CH;
                    Node.IsClusterHeads(j)=NON_CH;
                end
            end
        end
    EnRecP=EnRec(Elec,Kbit);
    Node.EnNode(j)=Node.EnNode(j)-EnRecP;
    Node.EnCCN(l,1)=Node.EnCCN(l,1)+EnRecP;
    end
end
CH=[];
CSET=zeros(NodeNums,NodeNums);%记录簇集成员节点ID
RCSET=zeros(NodeNums,NodeNums);
CHcount=length(find(Node.IsClusterHeads==2));
CH=find(Node.IsClusterHeads==2);
LBRU=zeros(1,NodeNums);
%%%%%%%%%%统计簇头竞争半径内的节点
for j=1:NodeNums
    if (Node.IsClusterHeads(j) == FINAL_CH)
        for z=1:NodeNums
            Djz=sqrt((Node.x(j)-Node.x(z))*(Node.x(j)-Node.x(z)) + (Node.y(j)-Node.y(z))*(Node.y(j)-Node.y(z)));
            if(Node.RC(j)>=Djz && j~=z)
               RCSET(j,z)=z;  %%保存簇头RC内的节点ID
            end
        end
    end
end

%%%%%%%%建立连接
for i=1:NodeNums
    if (sqrt(Node.x(i)^2 + Node.y(i)^2)<=do)
        Node.Cid(i)=0;%与基站连接
    else
        if(Node.IsClusterHeads(i)==NON_CH && Node.EnNode(i)>0)
           min_dis=Inf;%iInf表示正无穷大
           dist=Node.x(i)*Node.x(i)+Node.y(i)*Node.y(i);
           EntranPCH=EnTran(Elec,Eamp,Kbit,dist);
           Node.EnNode(i)=Node.EnNode(i)-EntranPCH;
           for c=1:1:CHcount
               temp=min(min_dis,sqrt((Node.x(CH(c))-Node.x(i))^2 + (Node.y(CH(c))-Node.y(i))^2));
               if (temp<min_dis)
                   min_dis=temp;
                   Node.Cid(i)=CH(c);%保存普通所连接的簇头的编号
               end
           end
           CSET(Node.Cid(i),i)=i;%保存对应簇头簇集内成员节点
           LBRU(Node.Cid(i))=LBRU(Node.Cid(i))+1;
        end
    end
end
%%%%%%%%%%%簇头在RC中选择最优传输路径
for i=1:NodeNums
    if(Node.IsClusterHeads(i)==FINAL_CH && Node.EnNode(i)>0)
       min_dis=Inf;%iInf表示正无穷大
       for c=1:1:length(Rc)
           temp=min(min_dis,(sqrt((Node.x(Rc(c))-Node.x(i))^2 + (Node.y(Rc(c))-Node.y(i))^2)+Node.x(Rc(c))^2+Node.y(Rc(c))^2));
           if (temp<min_dis)
               min_dis=temp;
               Node.Cid(i)=Rc(c);%保存普通所连接的簇头的编号
           end
       end
       LBRU(Node.Cid(i))=LBRU(Node.Cid(i))+1;
    end
end
%%%%%%%%计算普通节点加入簇集时的簇头接收能耗
for i=1:NodeNums
    Kbit=200;
    if(Node.IsClusterHeads(i) == FINAL_CH)
        Kbit=Kbit*LBRU(i);
        EnRecP=EnRec(Elec,Kbit);
        Node.EnNode(i)=Node.EnNode(i)-EnRecP;
        Node.EnCCN(l,1)=Node.EnCCN(l,1)+EnRecP;
    end
end

for r=1:1:rmax
dead=0;
r
for i=1:1:NodeNums
    if (Node.EnNode(i)<=0)
        dead=dead+1; 
        if(dead==1)
            if(Flag_FND(1)==0)
                FND(l,1)=r;
                Flag_FND(1)=1;
            end
        end
        if(dead==0.5*NodeNums)
            if(Flag_HND(1)==0)
                HND(l,1)=r;
                Flag_HND(1)=1;
            end
        end
    end
end
STATISTICS.DEAD(r)=dead;%死亡节点总数
STATISTICS.ALLIVE(r)=allive-dead;%存活节点数
%%%更新邻居节点
 for i=1:NodeNums
     if(Node.EnNode(i)>0)
         count =0 ;
        for j=1:NodeNums  
            if(j~=i && Node.EnNode(j)>0) 
            dist = ((Node.x(i)-Node.x(j)).^2)+((Node.y(i)-Node.y(j)).^2);  % the distance.^2
                   if dist < NodeTranR^2           %       考虑cost设置不同的时候需要将dist的值写入进去
                       count=count+1;
                       Node.Nbr(i,count)=j;       
                   end    
             end
             if j== NodeNums 
                    Node.NumNbr(i) = count ;          %邻居节点总个数
             end  
        end 
     end
 end                                               %初始化时求出节点的邻居节点及邻节点个数
if(mod(r,MaxInteral)==0) %%%轮换阶段
    ESET=zeros(NodeNums,NodeNums);%保存节点能量，便于在簇内对比
    for i=1:NodeNums
        if(Node.IsClusterHeads(i)==FINAL_CH)
            for j=1:NodeNums
                if(RCSET(i,j)>0 && Node.EnNode(RCSET(i,j))>0)
                    ESET(i,j)=Node.EnNode(RCSET(i,j));
                end
            end
        end
        Node.IsClusterHeads(i)=NON_CH;%重置簇头类型
    end
    [E,I]=max(ESET,[],2);
    
for i=1:length(I) %更新簇头
    Node.IsClusterHeads(I(i))=FINAL_CH;
    dist=Node.x(I(i))^2+Node.y(I(i)).^2;
    EntranPCH=EnTran(Elec,Eamp,Kbit,dist);
    Node.EnNode(i)=Node.EnNode(i)-EntranPCH;
    Node.EnCCN(l,1)=Node.EnCCN(l,1)+EntranPCH;
end
    
    
CH=[];
CSET=zeros(NodeNums,NodeNums);%记录簇集成员节点ID
CHcount=length(find(Node.IsClusterHeads==2));
CH=find(Node.IsClusterHeads==2);
LBRU=zeros(1,NodeNums);
Rc=[];
Rc=find(Node.IsClusterHeads==FINAL_CH);
%%%%%%%%重新建立连接
for i=1:NodeNums
    if (sqrt(Node.x(i)^2 + Node.y(i)^2)<=do)
        Node.Cid(i)=0;%与基站连接
    else
        if(Node.IsClusterHeads(i)==NON_CH && Node.EnNode(i)>0)
           min_dis=Inf;%iInf表示正无穷大
           dist=Node.x(i)*Node.x(i)+Node.y(i)*Node.y(i);
           EntranPCH=EnTran(Elec,Eamp,Kbit,dist);
           Node.EnNode(i)=Node.EnNode(i)-EntranPCH;
           for c=1:1:CHcount
               temp=min(min_dis,sqrt((Node.x(CH(c))-Node.x(i))^2 + (Node.y(CH(c))-Node.y(i))^2));
               if (temp<min_dis)
                   min_dis=temp;
                   Node.Cid(i)=CH(c);%保存普通所连接的簇头的编号
               end
           end   
           CSET(Node.Cid(i),i)=i;%保存对应簇头簇集内成员节点
           LBRU(Node.Cid(i))=LBRU(Node.Cid(i))+1;%计算簇内节点个数
        end
    end
end
%%%%%%%%%%%簇头在RC中选择最优传输路径
for i=1:NodeNums
    if(Node.IsClusterHeads(i)==FINAL_CH && Node.EnNode(i)>0)
       min_dis=Inf;%iInf表示正无穷大
       for c=1:1:length(Rc)
           temp=min(min_dis,(sqrt((Node.x(Rc(c))-Node.x(i))^2 + (Node.y(Rc(c))-Node.y(i))^2)+Node.x(Rc(c))^2+Node.y(Rc(c))^2));
           if (temp<min_dis)
               min_dis=temp;
               Node.Cid(i)=Rc(c);%保存普通所连接的簇头的编号
           end
       end
       LBRU(Node.Cid(i))=LBRU(Node.Cid(i))+1;
    end
end
%%%%%%%%计算普通节点加入簇集时的簇头接收能耗
for i=1:NodeNums
    Kbit=200;
    if(Node.IsClusterHeads(i) == FINAL_CH)
        Kbit=Kbit*LBRU(i);
        EnRecP=EnRec(Elec,Kbit);
        Node.EnNode(i)=Node.EnNode(i)-EnRecP;
        Node.EnCCN(l,1)=Node.EnCCN(l,1)+EnRecP;
    end
end
end
 
if(dead>0)
m=STATISTICS.DEAD(r-1);
if((dead-m)>0) %出现新死亡节点，在簇头集群内重新进行簇头筛选   
while sum(Node.Isstop)~=0
    iteration=iteration+1;
    for i =1:NodeNums             %此时进行选簇迭代过程，由于仿真每次都有节点id为1到NodeNums进行，导致id号靠前的节点其当选簇头概率较大，能耗较多。改进时可将节点id每次随机化
        if Node.Isstop(i)==1                     
          if Node.CHprob(i)<1         
             if Node.tent_CH(i)==NON_CH
              if rand(1,1)<Node.CHprob(i)
                 Node.IsClusterHeads(i)=TENTATIVE_CH; 
                 Node.tent_CH(i)=TOS_LOCAL_ADDRESS;
                 Node.tent_CH_Cost(i)=Node.NumNbr(i);                                  %cost值设置为节点的邻节点数目？即以度作为其归蔟度量值
              end
              %elseif  Node.tent_CH(i)==TOS_LOCAL_ADDRESS
             end
           Node.CHprob(i)=Node.CHprob(i).*2;
          else            %即CHprob等于1时的情况
             for j=1:NodeNums   %Node.n_finalCH(i)
                if Node.ListfinalCH(i,j) ~=0
                 if Node.my_final_CH_Cost(i) > Node.ListfinalCH_Cost(i,j) 
                     Node.my_finalCH(i)= Node.ListfinalCH(i,j);                           %进行归蔟行动
                     Node.my_final_CH_Cost(i)=Node.ListfinalCH_Cost(i,j);
                 end
                end  
             end
              % choose cluster head
            Node.Isstop(i)=0;   % diedai end   until CHprob==1
            if Node.my_finalCH(i) ~= NON_CH
                Node.IsClusterHeads(i)= NON_CH;                                                                      %这个为什么？
                Node.c(i)=Node.my_finalCH(i);
                Node.chcost=Node.my_final_CH_Cost(i);
                %join the cluster 
                dist =((Node.x(i)-Node.x(curentnbr)).^2)+((Node.y(i)-Node.y(curentnbr)).^2); % ((Node.x(i)-Node.x(curentnbr)).^2)+((Node.y(i)-Node.y(curentnbr)).^2);  % the distance.^2                这个值大概用来广播信息时用的
                EntranPCH=EnTran(Elec,Eamp,Kbit,dist) ;                                                          %大概是传播能量损耗值计算 函数
                Node.d(i)=((Node.x(i)-Node.x(Node.c(i))).^2)+((Node.y(i)-Node.y(Node.c(i))).^2);  % the distance.^2
                Node.EnNode(i)=Node.EnNode(i)-EntranPCH;                 %预先判断其能量是否足以进行这次传输
                if Node.EnNode(i) <= 0
                    Node.StateNode(i)=0;
                    Node.Isstop(i)=0;
                    Node.EnNode(i)=0;
                end
                EnRecP=EnRec(Elec,Kbit);                                          %这个是簇头接收数据后的能耗计算函数
                Node.EnNode(Node.c(i))=Node.EnNode(Node.c(i))-EnRecP;                
                if Node.EnNode(Node.c(i)) <= 0
                    Node.StateNode(Node.c(i))=0;
                    Node.Isstop(Node.c(i))=0;
                    Node.EnNode(Node.c(i))=0; 
                else                   
                    Node.csize(Node.c(i))=Node.csize(Node.c(i))+1;  % cluster size add one
                end
            else
                Node.IsClusterHeads(i)= FINAL_CH;
                Node.my_finalCH(i)=TOS_LOCAL_ADDRESS;
                Node.c(i)=TOS_LOCAL_ADDRESS;
                Node.my_final_CH_Cost(i)= Node.NumNbr(i);%computeDegree(i);
                Node.chcost=Node.my_final_CH_Cost(i);
                Node.d(i)=((Node.x(i)-Bx).^2)+((Node.y(i)-By).^2);  % the distance.^2  %簇头节点到sink节点的距离
                ClusterHeadNum=ClusterHeadNum+1;
            end
          end    
       end
      % compute consume energy
    if Node.IsClusterHeads(i) == TENTATIVE_CH  % & Node.tent_CH(curentnbr)==TOS_LOCAL_ADDRESS
        dist =NodeTranR.^2; % ((Node.x(i)-Node.x(curentnbr)).^2)+((Node.y(i)-Node.y(curentnbr)).^2);  % the distance.^2                    大概需要改正
           EntranPCH=EnTran(Elec,Eamp,Kbit,dist) ;                                 %即减去广播信息能耗
           Node.EnNode(i)=Node.EnNode(i)-EntranPCH;
           if Node.EnNode(i) <= 0
                Node.StateNode(i)=0;
                Node.Isstop(i)=0;
                Node.EnNode(i) =0;
           end
        for j=1:Node.NumNbr(i)
            curentnbr = Node.Nbr(i,j);                %处在临时簇头状态的蔟发送广播信息 求其邻居节点
             EnRecP=EnRec(Elec,Kbit);
           Node.EnNode(curentnbr)=Node.EnNode(curentnbr)-EnRecP;                  %这点是否需要加上在临时簇头节点能量值大于0的情况下
           if Node.EnNode(curentnbr) > 0
              if (Node.ListtentCH(curentnbr,i)==0)
                   Node.n_tentCH(curentnbr)=Node.n_tentCH(curentnbr)+1;    
              end 
              Node.ListtentCH(curentnbr,i) = i;
              Node.ListtentCH_Cost(curentnbr,i)=Node.NumNbr(i);%Node.computeDegree(i);

              % if Node.tent_CH(curentnbr)~=TOS_LOCAL_ADDRESS &(Node.tent_CH(curentnbr)==NON_CH | Node.tent_CH_Cost(i)< Node.tent_CH_Cost(curentnbr)  | ((Node.tent_CH_Cost(i)== Node.tent_CH_Cost(curentnbr) ) & i < Node.tent_CH(curentnbr)))
                if(Node.tent_CH(curentnbr)==NON_CH || Node.tent_CH_Cost(i)< Node.tent_CH_Cost(curentnbr)  || ((Node.tent_CH_Cost(i)== Node.tent_CH_Cost(curentnbr) ) && i < Node.tent_CH(curentnbr)))
                    Node.tent_CH_Cost(curentnbr)=Node.tent_CH_Cost(i);
                    Node.tent_CH(curentnbr)=i;
                end     % if(Node.tent_CH(curentnbr)==NON_CH | Node.tent_CH_Cost(i)< Node.tent_CH_Cost(curentnbr)
           else
                 Node.StateNode(curentnbr)=0;
                 Node.Isstop(curentnbr)=0;
                 Node.EnNode(curentnbr)=0; 
           end  
       end
            elseif  Node.IsClusterHeads(i) == FINAL_CH
                dist = NodeTranR.^2; %((Node.x(i)-Node.x(curentnbr)).^2)+((Node.y(i)-Node.y(curentnbr)).^2);  % the distance.^2
                   EntranPCH=EnTran(Elec,Eamp,Kbit,dist);
                   Node.EnNode(i)=Node.EnNode(i)-EntranPCH;
                   if Node.EnNode(i) <= 0
                        Node.StateNode(i)=0;
                        Node.Isstop(i)=0;
                        Node.EnNode(i)=0;
                   end  
            for j=1:Node.NumNbr(i)
                curentnbr = Node.Nbr(i,j);
                EnRecP=EnRec(Elec,Kbit);
                Node.EnNode(curentnbr)=Node.EnNode(curentnbr)-EnRecP;
               if Node.EnNode(curentnbr) > 0
                   if (Node.ListfinalCH(curentnbr,i)==0)                                   %这个是否跟上面一样错误
                       Node.n_finalCH(curentnbr)=Node.n_finalCH(curentnbr)+1;    
                   end
                  Node.ListfinalCH(curentnbr,i)=i;
                  Node.ListfinalCH_Cost(curentnbr,i)=Node.NumNbr(i);%Node.computeDegree(i);
               else
                 Node.EnNode(curentnbr)=0; 
                 Node.StateNode(curentnbr)=0;
                 Node.Isstop(curentnbr)=0;
               end    
            end
    end
     end  
end  
    %%%%%%%%%聚类
DRC=[];
Rc=[];
Rc=find(Node.IsClusterHeads==FINAL_CH);
for i=1:NodeNums
    %%%%%计算簇头与基站的最大值与最小值
    for j=1:NodeNums
        if Node.IsClusterHeads(i) == FINAL_CH
            DRC(j)=sqrt(Node.x(j)*Node.x(j)+Node.y(j)*Node.y(j));
        end
    end
end
%%%%%%计算竞争半径
for j=1:NodeNums
    if (Node.IsClusterHeads(j) == FINAL_CH)
        Node.RC(j)=(1-(max(DRC)-sqrt(Node.x(j)*Node.x(j)+Node.y(j)*Node.y(j)))/(3*(max(DRC)-min(DRC))))*RC0;
        dist=Node.x(j)*Node.x(j)+Node.y(j)*Node.y(j);
        EntranPCH=EnTran(Elec,Eamp,Kbit,dist);
        Node.EnNode(j)=Node.EnNode(j)-EntranPCH;
        Node.EnCCN(l,1)=Node.EnCCN(l,1)+EntranPCH;
    end
end
%%%%%利用竞争半径减少簇头节点数量
for j=1:NodeNums
    if Node.IsClusterHeads(j) == FINAL_CH
        for z=1:1:NodeNums
            if(Node.IsClusterHeads(z) == FINAL_CH && j~=z)
                df=sqrt((Node.x(j)-Node.x(z))*(Node.x(j)-Node.x(z)) + (Node.y(j)-Node.y(z))*(Node.y(j)-Node.y(z)));
                if (df>Node.RC(j) && df>Node.RC(z))
                    Node.IsClusterHeads(j)=FINAL_CH;
                    Node.IsClusterHeads(z)=FINAL_CH;
                else
                    Node.IsClusterHeads(z)=NON_CH;
                    Node.IsClusterHeads(j)=NON_CH;
                end
            end
        end
    EnRecP=EnRec(Elec,Kbit);
    Node.EnNode(j)=Node.EnNode(j)-EnRecP;
    Node.EnCCN(l,1)=Node.EnCCN(l,1)+EnRecP;
    end
end

CH=[];
CSET=zeros(NodeNums,NodeNums);%记录簇集成员节点ID
CHcount=length(find(Node.IsClusterHeads==2));
CH=find(Node.IsClusterHeads==2);
LBRU=zeros(1,NodeNums);
%%%%%%%%建立连接
for i=1:NodeNums
    if (sqrt(Node.x(i)^2 + Node.y(i)^2)<=do)
        Node.Cid(i)=0;%与基站连接
    else
        if(Node.IsClusterHeads(i)==NON_CH && Node.EnNode(i)>0)
           min_dis=Inf;%iInf表示正无穷大
           dist=Node.x(i)*Node.x(i)+Node.y(i)*Node.y(i);
           EntranPCH=EnTran(Elec,Eamp,Kbit,dist);
           Node.EnNode(i)=Node.EnNode(i)-EntranPCH;
           for c=1:1:CHcount
               temp=min(min_dis,sqrt((Node.x(CH(c))-Node.x(i))^2 + (Node.y(CH(c))-Node.y(i))^2));
               if (temp<min_dis)
                   min_dis=temp;
                   Node.Cid(i)=CH(c);%保存普通所连接的簇头的编号
                   LBRU(Node.Cid(i))=LBRU(Node.Cid(i))+1;%计算簇内节点个数
               end
           end   
           CSET(Node.Cid(i),i)=i;%保存对应簇头簇集内成员节点
        end
    end
end
%%%%%%%%%%%簇头在RC中选择最优传输路径
for i=1:NodeNums
    if(Node.IsClusterHeads(i)==FINAL_CH && Node.EnNode(i)>0)
       min_dis=Inf;%iInf表示正无穷大
       for c=1:1:length(Rc)
           temp=min(min_dis,(sqrt((Node.x(Rc(c))-Node.x(i))^2 + (Node.y(Rc(c))-Node.y(i))^2)+Node.x(Rc(c))^2+Node.y(Rc(c))^2));
           if (temp<min_dis)
               min_dis=temp;
               Node.Cid(i)=Rc(c);%保存普通所连接的簇头的编号
               LBRU(Node.Cid(i))=LBRU(Node.Cid(i))+1;
           end
       end
    end
end
%%%%%%%%计算普通节点加入簇集时的簇头接收能耗
for i=1:NodeNums
    Kbit=200;
    if(Node.IsClusterHeads(i) == FINAL_CH)
        Kbit=Kbit*LBRU(i);
        EnRecP=EnRec(Elec,Kbit);
        Node.EnNode(i)=Node.EnNode(i)-EnRecP;
        Node.EnCCN(l,1)=Node.EnCCN(l,1)+EnRecP;
    end
end
end
end
%%%%%%%%%%%%%%%%%%%%%能耗模型%%%%%%
for i=1:NodeNums
     if(Node.EnNode(i)>0)
         if(Node.Cid(i)~=0 )
             d=sqrt((Node.x(i)-Node.x(Node.Cid(i)))^2 + (Node.y(i)-Node.y(Node.Cid(i))^2));
            if(d>do)
               Node.EnNode(i)=Node.EnNode(i)-((EPC+ETX)*(4000)+(ERX+EDA)*(4000)*LBRU(i)+Emp*4000*(d * d * d * d));
            end
            if (d<=do)
                Node.EnNode(i)=Node.EnNode(i)-((EPC+ETX)*(4000)+(ERX+EDA)*(4000)*LBRU(i)+Efs*4000*(d * d)); 
            end
         else
            d=sqrt((Node.x(i))^2 + (Node.y(i))^2);
            if(d>do)
               Node.EnNode(i)=Node.EnNode(i)-((EPC+ETX)*(4000)+(ERX+EDA)*(4000)*LBRU(i)+Emp*4000*(d * d * d * d));
            end
            if (d<=do)
                Node.EnNode(i)=Node.EnNode(i)-((EPC+ETX)*(4000)+(ERX+EDA)*(4000)*LBRU(i)+Efs*4000*(d * d)); 
            end
         end
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%计算簇头能耗%%%%%%%%%%%%             
         if(Node.IsClusterHeads(i)==FINAL_CH)
             if(Node.Cid(i)~=0)
                d=sqrt((Node.x(i)-Node.x(Node.Cid(i)))^2 + (Node.y(i)-Node.y(Node.Cid(i))^2));
                if(d>do)
                   Node.EnCCN(l,1)=Node.EnCCN(l,1)+((EPC+ETX)*(4000)+(ERX+EDA)*(4000)*LBRU(i)+Emp*4000*(d * d * d * d));
                end
                if (d<=do)
                    Node.EnCCN(l,1)=Node.EnCCN(l,1)+((EPC+ETX)*(4000)+(ERX+EDA)*(4000)*LBRU(i)+Efs*4000*(d * d)); 
                end
             else
                d=sqrt((Node.x(i))^2 + (Node.y(i))^2);
                if(d>do)
                   Node.EnCCN(l,1)=Node.EnCCN(l,1)+((EPC+ETX)*(4000)+(ERX+EDA)*(4000)*LBRU(i)+Emp*4000*(d * d * d * d));
                end
                if (d<=do)
                    Node.EnCCN(l,1)=Node.EnCCN(l,1)+((EPC+ETX)*(4000)+(ERX+EDA)*(4000)*LBRU(i)+Efs*4000*(d * d)); 
                end
             end
         end
         
     end
end
Es=0;
for i=1:NodeNums
    if(Node.EnNode(i)>0)
        Es=Node.EnNode(i)+Es;
    end
end
Node.EnCCN(l,2)=Es;
for i=1:1:n
  if(Node.EnNode(i)<=0 && Node.AvL(l,i)==0)
    Node.AvL(l,i)=r;
    if(dead>=0.5*NodeNums)
        for j=1:NodeNums
            if(Node.AvL(l,j)==0)
                Node.AvL(l,j)=r;
            end
        end
    end
  end
end
for i=1:n
    if(Node.EnNode(i)>0 && dead<=0.5*n)
        data=data+1;
    end
end
end
Node.Data(1,l)=data;

%%%%%%%%%%%%%%%%%%%RE-HEED
%求邻节点个数
 for i=1:NodeNums
     count =0 ;
    for j=1:NodeNums  
        if(j~=i) 
        dist = ((Node2.x(i)-Node2.x(j)).^2)+((Node2.y(i)-Node2.y(j)).^2);  % the distance.^2
               if dist < NodeTranR^2           %       考虑cost设置不同的时候需要将dist的值写入进去
                   count=count+1;
                   Node2.Nbr(i,count)=j;       
               end    
         end
         if j== NodeNums 
                Node2.NumNbr(i) = count ;          %邻居节点总个数
         end  
    end 
 end                                               %初始化时求出节点的邻居节点及邻节点个数
 sym iteration;
for i=1:NodeNums
    Node2.CHprob(i)=Cprob*((Node2.EnNode(i))/MaxEn2(i));
end
ClusterHeadNum2=0;
iteration=0;
Node2.RC=zeros(1,NodeNums); 
dead2=0;
allive2=NodeNums;
%选簇
while sum(Node2.Isstop)~=0
    iteration=iteration+1;
    for i =1:NodeNums             %此时进行选簇迭代过程，由于仿真每次都有节点id为1到NodeNums进行，导致id号靠前的节点其当选簇头概率较大，能耗较多。改进时可将节点id每次随机化
        if Node2.Isstop(i)==1                      
          if Node2.CHprob(i)<1         
             if Node2.tent_CH(i)==NON_CH
              if rand(1,1)<Node2.CHprob(i)
                 Node2.IsClusterHeads(i)=TENTATIVE_CH; 
                 Node2.tent_CH(i)=TOS_LOCAL_ADDRESS;
                 Node2.tent_CH_Cost(i)=Node2.NumNbr(i);                                  %cost值设置为节点的邻节点数目？即以度作为其归蔟度量值
              end
              %elseif  Node.tent_CH(i)==TOS_LOCAL_ADDRESS
             end
           Node2.CHprob(i)=Node2.CHprob(i).*2;
          else            %即CHprob等于1时的情况
             for j=1:NodeNums   %Node.n_finalCH(i)
                if Node2.ListfinalCH(i,j) ~=0
                 if Node2.my_final_CH_Cost(i) > Node2.ListfinalCH_Cost(i,j) 
                     Node2.my_finalCH(i)= Node2.ListfinalCH(i,j);                           %进行归蔟行动
                     Node2.my_final_CH_Cost(i)=Node2.ListfinalCH_Cost(i,j);
                 end
                end  
             end
              % choose cluster head
              Node2.Isstop(i)=0;   % diedai end   until CHprob==1
              if Node2.my_finalCH(i) ~= NON_CH
                 Node2.IsClusterHeads(i)= NON_CH;                                                                      %这个为什么？
                 Node2.c(i)=Node2.my_finalCH(i);
                 Node2.chcost=Node2.my_final_CH_Cost(i);
                 %join the cluster 
                 dist =((Node2.x(i)-Node2.x(curentnbr)).^2)+((Node2.y(i)-Node2.y(curentnbr)).^2); 
                 EntranPCH=EnTran(Elec,Eamp,Kbit,dist);                                                          %大概是传播能量损耗值计算 函数
                 Node2.d(i)=((Node2.x(i)-Node2.x(Node2.c(i))).^2)+((Node2.y(i)-Node2.y(Node2.c(i))).^2);  % the distance.^2
                 Node2.EnNode(i)=Node2.EnNode(i)-EntranPCH;                 %预先判断其能量是否足以进行这次传输
                 if Node2.EnNode(i) <= 0
                    Node2.StateNode(i)=0;
                    Node2.Isstop(i)=0;
                    Node2.EnNode(i)=0;
                 end
                EnRecP=EnRec(Elec,Kbit);                                          %这个是簇头接收数据后的能耗计算函数
                Node2.EnNode(Node2.c(i))=Node2.EnNode(Node2.c(i))-EnRecP;                
                if Node2.EnNode(Node2.c(i)) <= 0
                    Node2.StateNode(Node2.c(i))=0;
                    Node2.Isstop(Node2.c(i))=0;
                    Node2.EnNode(Node2.c(i))=0; 
                else                   
                    Node2.csize(Node2.c(i))=Node2.csize(Node2.c(i))+1;  % cluster size add one
                end
             else
                 Node2.IsClusterHeads(i)= FINAL_CH;
                 Node2.my_finalCH(i)=TOS_LOCAL_ADDRESS;
                 Node2.c(i)=TOS_LOCAL_ADDRESS;
                 Node2.my_final_CH_Cost(i)= Node2.NumNbr(i);%computeDegree(i);
                 Node2.chcost=Node2.my_final_CH_Cost(i);
                 Node2.d(i)=((Node2.x(i)-Bx).^2)+((Node2.y(i)-By).^2);  % the distance.^2  %簇头节点到sink节点的距离
                 ClusterHeadNum2=ClusterHeadNum2+1;
              end
          end    
       end
      % compute consume energy
    if Node2.IsClusterHeads(i) == TENTATIVE_CH  % & Node.tent_CH(curentnbr)==TOS_LOCAL_ADDRESS
        dist =NodeTranR.^2; % ((Node.x(i)-Node.x(curentnbr)).^2)+((Node.y(i)-Node.y(curentnbr)).^2);  % the distance.^2                    大概需要改正
           EntranPCH=EnTran(Elec,Eamp,Kbit,dist) ;                                 %即减去广播信息能耗
           Node2.EnNode(i)=Node2.EnNode(i)-EntranPCH;
           if Node2.EnNode(i) <= 0
                Node2.StateNode(i)=0;
                Node2.Isstop(i)=0;
                Node2.EnNode(i) =0;
           end
        for j=1:Node2.NumNbr(i)
            curentnbr = Node2.Nbr(i,j);                %处在临时簇头状态的蔟发送广播信息 求其邻居节点
             EnRecP=EnRec(Elec,Kbit);
           Node2.EnNode(curentnbr)=Node2.EnNode(curentnbr)-EnRecP;                  %这点是否需要加上在临时簇头节点能量值大于0的情况下
           if Node2.EnNode(curentnbr) > 0
              if (Node2.ListtentCH(curentnbr,i)==0)
                   Node2.n_tentCH(curentnbr)=Node2.n_tentCH(curentnbr)+1;    
              end 
              Node2.ListtentCH(curentnbr,i) = i;
              Node2.ListtentCH_Cost(curentnbr,i)=Node2.NumNbr(i);%Node.computeDegree(i);

              % if Node.tent_CH(curentnbr)~=TOS_LOCAL_ADDRESS &(Node.tent_CH(curentnbr)==NON_CH | Node.tent_CH_Cost(i)< Node.tent_CH_Cost(curentnbr)  | ((Node.tent_CH_Cost(i)== Node.tent_CH_Cost(curentnbr) ) & i < Node.tent_CH(curentnbr)))
                if(Node2.tent_CH(curentnbr)==NON_CH || Node2.tent_CH_Cost(i)< Node2.tent_CH_Cost(curentnbr)  || ((Node2.tent_CH_Cost(i)== Node2.tent_CH_Cost(curentnbr) ) && i < Node2.tent_CH(curentnbr)))
                    Node2.tent_CH_Cost(curentnbr)=Node2.tent_CH_Cost(i);
                    Node2.tent_CH(curentnbr)=i;
                end     % if(Node.tent_CH(curentnbr)==NON_CH | Node.tent_CH_Cost(i)< Node.tent_CH_Cost(curentnbr)
           else
                 Node2.StateNode(curentnbr)=0;
                 Node2.Isstop(curentnbr)=0;
                 Node2.EnNode(curentnbr)=0; 
           end  
           end
            elseif  Node2.IsClusterHeads(i) == FINAL_CH
                dist = NodeTranR.^2; %((Node.x(i)-Node.x(curentnbr)).^2)+((Node.y(i)-Node.y(curentnbr)).^2);  % the distance.^2
                   EntranPCH=EnTran(Elec,Eamp,Kbit,dist);
                   Node2.EnNode(i)=Node2.EnNode(i)-EntranPCH;
                   if Node2.EnNode(i) <= 0
                        Node2.StateNode(i)=0;
                        Node2.Isstop(i)=0;
                        Node2.EnNode(i)=0;
                   end  
            for j=1:Node2.NumNbr(i)
                curentnbr = Node2.Nbr(i,j);
                EnRecP=EnRec(Elec,Kbit);
                Node2.EnNode(curentnbr)=Node2.EnNode(curentnbr)-EnRecP;
               if Node2.EnNode(curentnbr) > 0
                   if (Node2.ListfinalCH(curentnbr,i)==0)                                   %这个是否跟上面一样错误
                       Node2.n_finalCH(curentnbr)=Node2.n_finalCH(curentnbr)+1;    
                   end
                  Node2.ListfinalCH(curentnbr,i)=i;
                  Node2.ListfinalCH_Cost(curentnbr,i)=Node2.NumNbr(i);%Node.computeDegree(i);
               else
                 Node2.EnNode(curentnbr)=0; 
                 Node2.StateNode(curentnbr)=0;
                 Node2.Isstop(curentnbr)=0;
               end    
            end
        end
      end  
 end  
%%%%%%%%%聚类
DRC2=[];
Rc2=[];
Rc2=find(Node2.IsClusterHeads==FINAL_CH);
CSET2=zeros(NodeNums,NodeNums);%记录簇头RC半径内成员节点ID
for i=1:NodeNums
    %%%%%计算簇头与基站的最大值与最小值
    for j=1:NodeNums
        if Node2.IsClusterHeads(i) == FINAL_CH
            DRC2(j)=sqrt(Node2.x(j)*Node2.x(j)+Node2.y(j)*Node2.y(j));
        end
    end
end
%%%%%%计算竞争半径
for j=1:NodeNums
    if (Node2.IsClusterHeads(j) == FINAL_CH)
        Node2.RC(j)=(1-(max(DRC2)-sqrt(Node2.x(j)*Node2.x(j)+Node2.y(j)*Node2.y(j)))/(3*(max(DRC2)-min(DRC2))))*RC0;
        dist=Node2.x(j)*Node2.x(j)+Node2.y(j)*Node2.y(j);
        EntranPCH=EnTran(Elec,Eamp,Kbit,dist);
        Node2.EnNode(j)=Node2.EnNode(j)-EntranPCH;
        Node2.EnCCN(l,1)=Node2.EnCCN(l,1)+EntranPCH;
    end
end
Node2.RC=zeros(1,NodeNums)+max(Node2.RC);
%%%%%%利用竞争半径减少簇头节点数量
for j=1:NodeNums
    if Node2.IsClusterHeads(j) == FINAL_CH
        for z=1:1:NodeNums
            if(Node2.IsClusterHeads(z) == FINAL_CH && j~=z)
                df=sqrt((Node2.x(j)-Node2.x(z))*(Node2.x(j)-Node2.x(z)) + (Node2.y(j)-Node2.y(z))*(Node2.y(j)-Node2.y(z)));
                if (df>Node2.RC(j) && df>Node2.RC(z))
                    Node2.IsClusterHeads(j)=FINAL_CH;
                    Node2.IsClusterHeads(z)=FINAL_CH;
                else
                    Node2.IsClusterHeads(z)=NON_CH;
                    Node2.IsClusterHeads(j)=NON_CH;
                end
            end
        end
    EnRecP=EnRec(Elec,Kbit);
    Node2.EnNode(j)=Node2.EnNode(j)-EnRecP;
    Node2.EnCCN(l,1)=Node2.EnCCN(l,1)+EnRecP;
    end
end
%%%%%%%%%%统计簇头竞争半径内的节点
for j=1:NodeNums
    if (Node2.IsClusterHeads(j) == FINAL_CH)
        for z=1:NodeNums
            Djz=sqrt((Node2.x(j)-Node2.x(z))*(Node2.x(j)-Node2.x(z)) + (Node2.y(j)-Node2.y(z))*(Node2.y(j)-Node2.y(z)));
            if(Node2.RC(j)>=Djz && j~=z)
               CSET2(j,z)=z;
            end
        end
    end
end

CH2=[];
CHcount2=length(find(Node2.IsClusterHeads==2));
CH2=find(Node2.IsClusterHeads==2);
LBRU2=zeros(1,NodeNums);
%%%%%%%%建立连接
for i=1:NodeNums
    if (sqrt(Node2.x(i)^2 + Node2.y(i)^2)<=43.8529)
        Node2.Cid(i)=0;%与基站连接
    else
        if(Node2.IsClusterHeads(i)==NON_CH && Node2.EnNode(i)>0)
           min_dis=Inf;%iInf表示正无穷大
           dist=Node2.x(i)*Node2.x(i)+Node2.y(i)*Node2.y(i);
           EntranPCH=EnTran(Elec,Eamp,Kbit,dist);
           Node2.EnNode(i)=Node2.EnNode(i)-EntranPCH;
           for c=1:1:CHcount2
               temp=min(min_dis,sqrt((Node2.x(CH2(c))-Node2.x(i))^2 + (Node2.y(CH2(c))-Node2.y(i))^2));
               if (temp<min_dis)
                   min_dis=temp;
                   Node2.Cid(i)=CH2(c);%保存普通所连接的簇头的编号
               end
           end
           LBRU2(Node2.Cid(i))=LBRU2(Node2.Cid(i))+1;
        end
    end
end
%%%%%%%%%%%簇头在RC中选择最优传输路径
for i=1:NodeNums
    if(Node2.IsClusterHeads(i)==FINAL_CH && Node2.EnNode(i)>0)
       min_dis=Inf;%iInf表示正无穷大
       for c=1:1:length(Rc2)
           temp=min(min_dis,(sqrt((Node2.x(Rc2(c))-Node2.x(i))^2 + (Node2.y(Rc2(c))-Node2.y(i))^2)+Node2.x(Rc2(c))^2+Node2.y(Rc2(c))^2));
           if (temp<min_dis)
               min_dis=temp;
               Node2.Cid(i)=Rc2(c);%保存普通所连接的簇头的编号
           end
       end
       LBRU2(Node2.Cid(i))=LBRU2(Node2.Cid(i))+1;
    end
end
%%%%%%%%计算普通节点加入簇集时的簇头接收能耗
for i=1:NodeNums
    Kbit=200;
    if(Node2.IsClusterHeads(i) == FINAL_CH)
        Kbit=Kbit*LBRU1(i);
        EnRecP=EnRec(Elec,Kbit);
        Node2.EnNode(i)=Node2.EnNode(i)-EnRecP;
        Node2.EnCCN(l,1)=Node2.EnCCN(l,1)+EnRecP;
    end
end

for r=1:1:rmax
dead2=0;
r
for i=1:1:NodeNums
    if (Node2.EnNode(i)<=0)
        dead2=dead2+1;
        if(dead2==1)
            if(Flag_FND(3)==0)
                FND(l,3)=r;
                Flag_FND(3)=1;
            end
        end
        if(dead2==0.5*NodeNums)
            if(Flag_HND(3)==0)
                HND(l,3)=r;
                Flag_HND(3)=1;
            end
        end
    end
end
STATISTICS.DEAD2(r)=dead2;%死亡节点总数
STATISTICS.ALLIVE2(r)=allive2-dead2;%存活节点数
m=allive2-dead2;
if(dead2>0)
m2=STATISTICS.DEAD(r-1);
if((dead2-m2)>0) %出现新死亡节点，在簇头竞争半径集群内重新进行簇头筛选
    ESET2=zeros(NodeNums,NodeNums);%保存节点能量，便于在簇内对比
    for i=1:NodeNums
        if(Node2.IsClusterHeads(i)==FINAL_CH)
            for j=1:NodeNums
                if(CSET2(i,j)>0 && Node2.EnNode(CSET2(i,j))>0)
                    ESET2(i,j)=Node2.EnNode(CSET2(i,j));
                end
            end
        end
        Node2.IsClusterHeads(i)=NON_CH;%重置簇头类型
    end
    [E,I]=max(ESET2,[],2);
    for i=1:length(I) %更新簇头
        Node2.IsClusterHeads(I(i))=FINAL_CH;
        dist=Node2.x(I(i))^2+Node2.y(I(i)).^2;
        EntranPCH=EnTran(Elec,Eamp,Kbit,dist);
        Node2.EnNode(i)=Node2.EnNode(i)-EntranPCH;
        Node2.EnCCN(l,1)=Node2.EnCCN(l,1)+EntranPCH;
    end
    %%%%%%%%%聚类
DRC2=[];
Rc2=[];
Rc2=find(Node2.IsClusterHeads==FINAL_CH);
for i=1:NodeNums
    %%%%%计算簇头与基站的最大值与最小值
    for j=1:NodeNums
        if Node2.IsClusterHeads(i) == FINAL_CH
            DRC2(j)=sqrt(Node2.x(j)*Node2.x(j)+Node2.y(j)*Node2.y(j));
        end
    end
    %%%%%%计算竞争半径
    for j=1:NodeNums
        if (Node2.IsClusterHeads(j) == FINAL_CH)
            Node2.RC(j)=(1-(max(DRC2)-sqrt(Node2.x(j)*Node2.x(j)+Node2.y(j)*Node2.y(j)))/(3*(max(DRC2)-min(DRC2))))*RC0;
            dist=Node2.x(j)*Node2.x(j)+Node2.y(j)*Node2.y(j);
            EntranPCH=EnTran(Elec,Eamp,Kbit,dist);
            Node2.EnNode(j)=Node2.EnNode(j)-EntranPCH;
            Node2.EnCCN(l,1)=Node2.EnCCN(l,1)+EntranPCH;
        end
    end
   Node2.RC=zeros(1,NodeNums)+median(Node2.RC);
    %%%%%%利用竞争半径减少簇头节点数量
    for j=1:NodeNums
        if Node2.IsClusterHeads(j) == FINAL_CH
            for z=1:1:NodeNums
                if(Node2.IsClusterHeads(z) == FINAL_CH && j~=z)
                    df=sqrt((Node2.x(j)-Node2.x(z))*(Node2.x(j)-Node2.x(z)) + (Node2.y(j)-Node2.y(z))*(Node2.y(j)-Node2.y(z)));
                    if (df>Node2.RC(j) && df>Node2.RC(z))
                        Node2.IsClusterHeads(j)=FINAL_CH;
                        Node2.IsClusterHeads(z)=FINAL_CH;
                    else
                        Node2.IsClusterHeads(z)=NON_CH;
                        Node2.IsClusterHeads(j)=NON_CH;
                    end
                end
            end
        EnRecP=EnRec(Elec,Kbit);
        Node2.EnNode(j)=Node2.EnNode(j)-EnRecP;
        Node2.EnCCN(l,1)=Node2.EnCCN(l,1)+EnRecP;
        end
    end 
    %%%%%%%%%%统计簇头竞争半径内的节点
    for j=1:NodeNums
        if (Node2.IsClusterHeads(j) == FINAL_CH)
            for z=1:NodeNums
                Djz=sqrt((Node2.x(j)-Node2.x(z))*(Node2.x(j)-Node2.x(z)) + (Node2.y(j)-Node2.y(z))*(Node2.y(j)-Node2.y(z)));
                if(Node2.RC(j)>=Djz && j~=z)
                   CSET2(j,z)=z;
                end
            end
        end
    end
end
CH2=[];
CSET2=zeros(NodeNums,NodeNums);%记录簇集成员节点ID
CHcount2=length(find(Node2.IsClusterHeads==2));
CH2=find(Node2.IsClusterHeads==2);
LBRU2=zeros(1,NodeNums);
%%%%%%%%建立连接
for i=1:NodeNums
    if (sqrt(Node2.x(i)^2 + Node2.y(i)^2)<=do)
        Node2.Cid(i)=0;%与基站连接
    else
        if(Node2.IsClusterHeads(i)==NON_CH && Node2.EnNode(i)>0)
           min_dis=Inf;%iInf表示正无穷大
           dist=Node2.x(i)*Node2.x(i)+Node2.y(i)*Node2.y(i);
           EntranPCH=EnTran(Elec,Eamp,Kbit,dist);
           Node2.EnNode(i)=Node2.EnNode(i)-EntranPCH;
           for c=1:1:CHcount2
               temp=min(min_dis,sqrt((Node2.x(CH2(c))-Node2.x(i))^2 + (Node2.y(CH2(c))-Node2.y(i))^2));
               if (temp<min_dis)
                   min_dis=temp;
                   Node2.Cid(i)=CH2(c);%保存普通所连接的簇头的编号
               end
           end
           CSET2(Node2.Cid(i),i)=i;%保存对应簇头簇集内成员节点
           LBRU2(Node2.Cid(i))=LBRU2(Node2.Cid(i))+1;%计算簇内节点个数
        end
    end
end
%%%%%%%%%%%簇头在RC中选择最优传输路径
for i=1:NodeNums
    if(Node2.IsClusterHeads(i)==FINAL_CH && Node2.EnNode(i)>0)
       min_dis=Inf;%iInf表示正无穷大
       for c=1:1:length(Rc2)
           temp=min(min_dis,(sqrt((Node2.x(Rc2(c))-Node2.x(i))^2 + (Node2.y(Rc2(c))-Node2.y(i))^2)+Node2.x(Rc2(c))^2+Node2.y(Rc2(c))^2));
           if (temp<min_dis)
               min_dis=temp;
               Node2.Cid(i)=Rc2(c);%保存普通所连接的簇头的编号
           end
       end
       LBRU2(Node2.Cid(i))=LBRU2(Node2.Cid(i))+1;
    end
end
end
end

%%%%%%%%计算普通节点加入簇集时的簇头接收能耗
for i=1:NodeNums
    Kbit=200;
    if(Node2.IsClusterHeads(i) == FINAL_CH)
        Kbit=Kbit*LBRU1(i);
        EnRecP=EnRec(Elec,Kbit);
        Node2.EnNode(i)=Node2.EnNode(i)-EnRecP;
        Node2.EnCCN(l,1)=Node2.EnCCN(l,1)+EnRecP;
    end
    eav=(ERX*4000*(sum(CSET2(i,:)>0,2)-1)+EDA*4000*sum(CSET2(i,:)>0,2))/(2*sum(CSET2(i,:)>0,2));
    for j=1:NodeNums
        if(CSET2(i,j)~=0)
            Node2.EnNode(CSET2(i,j))=Node2.EnNode(CSET2(i,j))-eav;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%能耗模型%%%%%%
for i=1:NodeNums
     if(Node2.EnNode(i)>0)
         if(Node2.IsClusterHeads(i) == FINAL_CH)
             d=sqrt((Node2.x(i)-Node2.x(Node2.Cid(i)))^2 + (Node2.y(i)-Node2.y(Node2.Cid(i))^2));
            if(d>do)
               Node2.EnNode(i)=Node2.EnNode(i)-((EPC+ETX)*(4000)+Emp*4000*(d * d * d * d));
               Node2.EnCCN(l,1)=Node2.EnCCN(l,1)+((EPC+ETX)*(4000)+Emp*4000*(d * d * d * d));
            end
            if (d<=do)
                Node2.EnNode(i)=Node2.EnNode(i)-((EPC+ETX)*(4000)+Efs*4000*(d * d));
                Node2.EnCCN(l,1)=Node2.EnCCN(l,1)+((EPC+ETX)*(4000)+Efs*4000*(d * d));
            end
         else
             if(Node2.Cid(i)~=0 )
                 d=sqrt((Node2.x(i)-Node2.x(Node2.Cid(i)))^2 + (Node2.y(i)-Node2.y(Node2.Cid(i))^2));
                if(d>do)
                   Node2.EnNode(i)=Node2.EnNode(i)-((EPC+ETX)*(4000)+(ERX+EDA)*(4000)*LBRU2(i)+Emp*4000*(d * d * d * d));
                end
                if (d<=do)
                    Node2.EnNode(i)=Node2.EnNode(i)-((EPC+ETX)*(4000)+(ERX+EDA)*(4000)*LBRU2(i)+Efs*4000*(d * d)); 
                end
             else
                d=sqrt((Node2.x(i))^2 + (Node2.y(i))^2);
                if(d>do)
                   Node2.EnNode(i)=Node2.EnNode(i)-((EPC+ETX)*(4000)+(ERX+EDA)*(4000)*LBRU2(i)+Emp*4000*(d * d * d * d));
                end
                if (d<=do)
                    Node2.EnNode(i)=Node2.EnNode(i)-((EPC+ETX)*(4000)+(ERX+EDA)*(4000)*LBRU2(i)+Efs*4000*(d * d)); 
                end
             end
         end
     end
end
Es=0;
for i=1:NodeNums
    if(Node2.EnNode(i)>0)
        Es=Node2.EnNode(i)+Es;
    end
end
Node2.EnCCN(l,2)=Es;
for i=1:1:n
  if(Node2.EnNode(i)<=0 && Node2.AvL(l,i)==0)
    Node2.AvL(l,i)=r;
    if(dead2>=0.5*NodeNums)
        for j=1:NodeNums
            if(Node2.AvL(l,j)==0)
                Node2.AvL(l,j)=r;
            end
        end
    end
  end
end
for i=1:n
    if(Node2.EnNode(i)>0 && dead2<=0.5*n)
        data2=data2+1;
    end
end
end
Node2.Data(1,l)=data2;
end

r=1:rmax;
STATISTICS.LIVE=allive-STATISTICS.DEAD;
STATISTICS.LIVE1=allive-STATISTICS.DEAD1;
STATISTICS.LIVE2=allive-STATISTICS.DEAD2;
STATISTICS.LIVE3=allive-STATISTICS.DEAD3;
figure('name','节点存活情况比较')
plot(r,STATISTICS.LIVE,'-b',r,STATISTICS.LIVE1,'-r',r,STATISTICS.LIVE2,'-k',r,STATISTICS.LIVE3,'-g');
legend('RUHEED','ASHEED','RE-HEED','LEACH-N');
xlabel('x(time)');
ylabel('y(live)');
title('\bfRUHEED、RUN、RE-HEED、LEACH-N存活节点数随时间的变化对比');

z=1:10;
figure('name','Total data transmission per unit of energy, NC=100')
plot(z,Node.Data,'-b',z,Node1.Data,'-r',z,Node2.Data,'-k',z,Node3.Data,'-g');
legend('RUHEED','ASHEED','ER-HEED','LEACH-N');
xlabel('X(Time)');
ylabel('Y(Data)');
title('\bfTotal data transfer per unit energy for RUHEED, ASHEED, ER-HEED, LEACH-N, NC=100');