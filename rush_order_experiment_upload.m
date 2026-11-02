clear; clc; rng(20260926,'twister');
assert(license('test','optimization_toolbox')==1,'Optimization Toolbox required.');
root=fileparts(mfilename('fullpath'));if isempty(root),root=pwd;end
out=fullfile(root,'rush_order_extended_results');if ~exist(out,'dir'),mkdir(out);end
runMode='paper'; 
cfg.pickers=3;cfg.orders=[6 8 10];
if strcmp(runMode,'paper')
 cfg.repetitions=15;cfg.maxSolverSeconds=30;
elseif strcmp(runMode,'pilot')
 cfg.repetitions=2;cfg.maxSolverSeconds=8;
else
 error('runMode must be paper or pilot.');
end
cfg.relativeGap=.05;
cfg.regularSlack=[-2 8]; 
cfg.rushSlack=[.05 2 8]; 
cfg.weights=[.25 .45 .20 .10]; 
cfg.reference=[35 120 100 3];cfg.maxReassigned=4;cfg.speed=60;

R=table();scenarioInputs=table();orderInputs=table();pickerInputs=table();
caseID=0;target=numel(cfg.orders)*cfg.repetitions* ...
 numel(cfg.regularSlack)*numel(cfg.rushSlack);
for ni=1:numel(cfg.orders)
 for rep=1:cfg.repetitions
  base=generateBase(cfg,cfg.orders(ni));baseCase=(ni-1)*cfg.repetitions+rep;
  for dl=1:numel(cfg.regularSlack)
   for ul=1:numel(cfg.rushSlack)
    caseID=caseID+1; I=base;
    I.due=max(I.t0+1,I.baselineCompletion(1:I.n)+ ...
      cfg.regularSlack(dl)+3*I.dueNoise);
    I.rushDue=I.earliestRush+cfg.rushSlack(ul);
    [si,oi,pi]=describeInputs(I,cfg,caseID,baseCase,ni,dl,ul,rep);
    scenarioInputs=[scenarioInputs;si];orderInputs=[orderInputs;oi]; 
    pickerInputs=[pickerInputs;pi];
    A=cell(3,1);A{1}=baselineAppend(I,cfg);
    A{2}=bestInsertion(I,cfg);A{3}=optimizeMILP(I,cfg,A{2});
    for strategyIndex=1:3
     a=A{strategyIndex};k=measure(I,a,cfg);
     row=table(caseID,baseCase,I.n,dl,ul,rep, ...
       string(['S' num2str(strategyIndex-1)]), ...
       k.response,k.tardiness,k.late,k.distance,k.added, ...
       k.reassigned,k.inversions,k.deadlineMiss,k.objective, ...
       a.runtime,a.gap,string(a.status), ...
       'VariableNames',{'instance','baseCase','regularOrders', ...
       'regularDeadlineLevel','rushDeadlineLevel','replicate','strategy', ...
       'rushResponse_min','regularTardiness_min','lateRegularOrders', ...
       'totalDistance_m','additionalDistance_m','reassignedOrders', ...
       'sequenceInversions','rushDeadlineMiss','objective','runtime_s', ...
       'relativeGap','solverStatus'});
     R=[R;row]; 
    end
   
    writetable(R,fullfile(out,'instance_results.csv'));
    writetable(scenarioInputs,fullfile(out,'scenario_inputs.csv'));
    writetable(orderInputs,fullfile(out,'order_inputs.csv'));
    writetable(pickerInputs,fullfile(out,'picker_inputs.csv'));
    fprintf('%d/%d | n=%d base=%d regular=%d rush=%d | S2 %s %.1f s\n', ...
      caseID,target,I.n,baseCase,dl,ul,A{3}.status,A{3}.runtime);
   end
  end
 end
end
metrics={'rushResponse_min','regularTardiness_min','lateRegularOrders', ...
 'totalDistance_m','additionalDistance_m','reassignedOrders', ...
 'sequenceInversions','rushDeadlineMiss','objective','runtime_s'};
summary=table();paired=table();pairedSummary=table();
for ni=1:numel(cfg.orders)
 for dl=1:numel(cfg.regularSlack)
  for ul=1:numel(cfg.rushSlack)
   subset=R(R.regularOrders==cfg.orders(ni) & ...
    R.regularDeadlineLevel==dl & R.rushDeadlineLevel==ul,:);
   for s=0:2
    t=subset(subset.strategy==string(['S' num2str(s)]),:);
    for q=1:numel(metrics)
     vals=t.(metrics{q});[lo,hi]=bootstrapMeanCI(vals,1000);
     row=table(cfg.orders(ni),dl,ul,string(['S' num2str(s)]), ...
       string(metrics{q}),numel(vals),mean(vals),std(vals),lo,hi, ...
       'VariableNames',{'regularOrders','regularDeadlineLevel', ...
       'rushDeadlineLevel','strategy','metric','independentBaseCases', ...
       'mean','sampleSD','bootstrap95Low','bootstrap95High'});
     summary=[summary;row]; 
    end
   end
   ids=unique(subset.baseCase);
   for b=1:2
    vals=zeros(numel(ids),3);
    for j=1:numel(ids)
     t=subset(subset.baseCase==ids(j),:);
     z=t(t.strategy=='S2',:);ref=t(t.strategy==string(['S' num2str(b-1)]),:);
     vals(j,:)=[z.objective-ref.objective, ...
       z.rushResponse_min-ref.rushResponse_min, ...
       z.regularTardiness_min-ref.regularTardiness_min];
     paired=[paired;table(ids(j),cfg.orders(ni),dl,ul, ...
       string(['S' num2str(b-1)]),vals(j,1),vals(j,2),vals(j,3), ...
       'VariableNames',{'baseCase','regularOrders','regularDeadlineLevel', ...
       'rushDeadlineLevel','comparator','deltaObjective', ...
       'deltaRush_min','deltaTardiness_min'})]; 
    end
    for q=1:3
     names={'deltaObjective','deltaRush_min','deltaTardiness_min'};
     [lo,hi]=bootstrapMeanCI(vals(:,q),1000);
     pairedSummary=[pairedSummary;table(cfg.orders(ni),dl,ul, ...
       string(['S' num2str(b-1)]),string(names{q}),numel(ids), ...
       mean(vals(:,q)),std(vals(:,q)),lo,hi, ...
       'VariableNames',{'regularOrders','regularDeadlineLevel', ...
       'rushDeadlineLevel','comparator','metric','independentBaseCases', ...
       'meanDifference','sampleSD','bootstrap95Low','bootstrap95High'})]; 
    end
   end
  end
 end
end
writetable(summary,fullfile(out,'scenario_summary.csv'));
writetable(paired,fullfile(out,'paired_differences.csv'));
writetable(pairedSummary,fullfile(out,'paired_summary.csv'));
factorEffects=table();
for ni=1:numel(cfg.orders)
 for s=0:2
  for dl=1:numel(cfg.regularSlack)
   for ul=2:numel(cfg.rushSlack)
    for metric=string({'rushResponse_min','regularTardiness_min','objective'})
     v=pairedFactorValues(R,cfg.orders(ni),s,dl,1,dl,ul,metric);
     [lo,hi]=bootstrapMeanCI(v,1000);
     factorEffects=[factorEffects;table(string('rushDeadline'), ...
       cfg.orders(ni),string(['S' num2str(s)]),dl,ul, ...
       metric,numel(v),mean(v),lo,hi, ...
       'VariableNames',{'factor','regularOrders','strategy', ...
       'fixedLevel','comparedLevel','metric','independentBaseCases', ...
       'meanDifference','bootstrap95Low','bootstrap95High'})]; 
    end
   end
  end
  for ul=1:numel(cfg.rushSlack)
   for metric=string({'rushResponse_min','regularTardiness_min','objective'})
    v=pairedFactorValues(R,cfg.orders(ni),s,1,ul,2,ul,metric);
    [lo,hi]=bootstrapMeanCI(v,1000);
    factorEffects=[factorEffects;table(string('regularDeadline'), ...
      cfg.orders(ni),string(['S' num2str(s)]),ul,2, ...
      metric,numel(v),mean(v),lo,hi, ...
      'VariableNames',{'factor','regularOrders','strategy', ...
      'fixedLevel','comparedLevel','metric','independentBaseCases', ...
      'meanDifference','bootstrap95Low','bootstrap95High'})]; 
   end
  end
 end
end
writetable(factorEffects,fullfile(out,'deadline_effects.csv'));
fig=figure('Visible','off','Color','w','Position',[100 100 1000 440]);
tiledlayout(1,2,'TileSpacing','compact');
nexttile;hold on;
for dl=1:numel(cfg.regularSlack)
 for ul=1:numel(cfg.rushSlack)
  t=pairedSummary(pairedSummary.regularDeadlineLevel==dl & ...
    pairedSummary.rushDeadlineLevel==ul & ...
    pairedSummary.comparator=='S1' & ...
    pairedSummary.metric=='deltaObjective',:);
  plot(t.regularOrders,t.meanDifference,'-o','DisplayName', ...
   sprintf('Regular %d, rush %d',dl,ul));
 end
end
xlabel('Regular orders');ylabel('S2 minus S1 objective');grid on;legend('Location','best');
nexttile;hold on;
for s=1:2
 for ul=1:numel(cfg.rushSlack)
  t=summary(summary.regularDeadlineLevel==1 & ...
    summary.rushDeadlineLevel==ul & ...
    summary.strategy==string(['S' num2str(s)]) & ...
    summary.metric=='rushResponse_min',:);
  plot(t.regularOrders,t.mean,'-o','DisplayName',sprintf('S%d, rush %d',s,ul));
 end
end
xlabel('Regular orders');ylabel('Rush response (min)');grid on;legend('Location','best');
exportgraphics(fig,fullfile(out,'scenario_comparison.png'),'Resolution',200);close(fig);
save(fullfile(out,'experiment.mat'),'R','summary','paired', ...
 'pairedSummary','factorEffects','scenarioInputs','orderInputs','pickerInputs','cfg');
fid=fopen(fullfile(out,'experiment_settings.txt'),'w');
fprintf(fid,['Mode: %s\nSeed: 20260926\nIndependent base cases per order count: %d\n', ...
 'Regular orders: %s; regular due slack (min): %s; rush due slack (min): %s\n', ...
 'Factorial design: %d independent base states x %d deadline conditions = %d instances.\n', ...
 '3 pickers; 100 x 100 m square; Euclidean travel; speed %.1f m/min; ', ...
 'service 2-6 min; regular due = max(1, baseline completion + regular slack ', ...
 '+ 3 x fixed per-order normal deviate).\n', ...
 'Weights: %s; fixed reference values: %s; Rmax: %d.\n', ...
 'S2 time limit: %.1f seconds; relative gap target: %.3f.\n', ...
 '95%% percentile bootstrap of means across independent base states, ', ...
 '1000 resamples. Interpret only per factorial cell; repeated variants ', ...
 'of one base state are dependent.\n', ...
 'Sequence inversions are a descriptive diagnostic, not in the objective.\n', ...
 'S0 may miss rush deadline; S1 and S2 enforce it. ', ...
 'S2 fallback retains S1 if the MILP incumbent is worse.\n'], ...
 runMode,cfg.repetitions,mat2str(cfg.orders),mat2str(cfg.regularSlack), ...
 mat2str(cfg.rushSlack),numel(cfg.orders)*cfg.repetitions, ...
 numel(cfg.regularSlack)*numel(cfg.rushSlack),target,cfg.speed, ...
 mat2str(cfg.weights),mat2str(cfg.reference),cfg.maxReassigned, ...
 cfg.maxSolverSeconds,cfg.relativeGap);
fclose(fid);disp('Finished. Scenario and paired summaries are ready.');

function [lo,hi]=bootstrapMeanCI(x,B)
 x=x(:);n=numel(x);assert(n>=2,'At least two independent cases needed.');
 idx=randi(n,n,B);samples=sort(mean(x(idx),1));
 lo=samples(max(1,ceil(.025*B)));
 hi=samples(min(B,ceil(.975*B)));
end
function v=pairedFactorValues(R,n,s,dl0,ul0,dl1,ul1,metric)
 t=R(R.regularOrders==n & R.strategy==string(['S' num2str(s)]),:);
 ids=unique(t.baseCase);v=zeros(numel(ids),1);
 for j=1:numel(ids)
  a=t(t.baseCase==ids(j) & t.regularDeadlineLevel==dl0 & ...
    t.rushDeadlineLevel==ul0,:);
  b=t(t.baseCase==ids(j) & t.regularDeadlineLevel==dl1 & ...
    t.rushDeadlineLevel==ul1,:);
  assert(height(a)==1 && height(b)==1,'Factorial pairing is incomplete.');
  v(j)=b.(char(metric))-a.(char(metric));
 end
end
function I=generateBase(c,n)
 I.n=n;I.P=c.pickers;I.start=100*rand(c.pickers,2);
 I.xy=100*rand(n+1,2);I.service=2+4*rand(n+1,1);
 I.t0=0;I.base=cell(c.pickers,1);I.owner=zeros(n,1);
 perm=randperm(n);
 for j=1:n
  p=mod(j-1,c.pickers)+1;I.base{p}(end+1)=perm(j);I.owner(perm(j))=p;
 end
 I.ds=zeros(c.pickers,n+1);I.d=zeros(n+1);
 for p=1:c.pickers,for j=1:n+1
  I.ds(p,j)=norm(I.start(p,:)-I.xy(j,:));end,end
 for i=1:n+1,for j=1:n+1
  I.d(i,j)=norm(I.xy(i,:)-I.xy(j,:));end,end
 [I.baselineCompletion,~,I.baselineDistance]=simulate(I,I.base,c);
 I.dueNoise=randn(n,1);
 I.earliestRush=min(I.ds(:,n+1)/c.speed)+I.service(n+1);
end
function [scenario,orders,pickers]=describeInputs(I,c,id,baseID,ni,dl,ul,rep)
 n=I.n;P=I.P;position=zeros(n,1);
 for p=1:P,position(I.base{p})=(1:numel(I.base{p}))';end
 scenario=table(id,baseID,ni,n,dl,ul,rep,c.regularSlack(dl), ...
   c.rushSlack(ul),I.t0,I.earliestRush,I.rushDue, ...
   I.baselineDistance,c.speed,P, ...
   'VariableNames',{'instance','baseCase','orderCountLevel', ...
    'regularOrders','regularDeadlineLevel','rushDeadlineLevel', ...
    'replicate','regularSlack_min','rushSlack_min','arrival_min', ...
    'earliestRushCompletion_min','rushDeadline_min', ...
    'baselineDistance_m','speed_m_per_min','pickers'});
 j=(1:n+1)';rush=j==n+1;due=[I.due;I.rushDue];
 orders=table(repmat(id,n+1,1),repmat(baseID,n+1,1),j,rush, ...
   I.xy(:,1),I.xy(:,2),I.service,due,[I.owner;NaN], ...
   [position;NaN], ...
   'VariableNames',{'instance','baseCase','orderID','isRush', ...
   'x_m','y_m','service_min','due_min','baselinePicker','baselineSequence'});
 pickers=table(repmat(id,P,1),repmat(baseID,P,1),(1:P)', ...
   I.start(:,1),I.start(:,2), ...
   'VariableNames',{'instance','baseCase','pickerID','startX_m','startY_m'});
end
function [C,D,dist,finish]=simulate(I,routes,c)
 C=zeros(I.n+1,1); D=zeros(I.P,1); finish=zeros(I.P,1);
 for p=1:I.P
  prev=0; t=I.t0;
  for j=routes{p}
   if prev==0, delta=I.ds(p,j); else, delta=I.d(prev,j); end
   D(p)=D(p)+delta; t=t+delta/c.speed+I.service(j);
   C(j)=t; prev=j;
  end
  finish(p)=t;
 end
 dist=sum(D);
end
function a=baselineAppend(I,c)
 tic; routes=I.base; [~,~,~,times]=simulate(I,routes,c);
 [~,p]=min(times); routes{p}(end+1)=I.n+1;
 a=pack(routes,toc,'rule',NaN);
end
function a=bestInsertion(I,c)
 tic; best=inf; routes=I.base;
 for p=1:I.P
  for pos=0:numel(I.base{p})
   r=I.base; seq=r{p}; r{p}=[seq(1:pos) I.n+1 seq(pos+1:end)];
   [completion,~,~]=simulate(I,r,c);
   if completion(end)>I.rushDue+1e-8, continue; end
   k=measure(I,pack(r,0,'enumerated',0),c);
   if k.objective<best, best=k.objective; routes=r; end
  end
 end
 if isinf(best)
  a=baselineAppend(I,c);a.status='no feasible insertion';a.runtime=toc;
 else
  a=pack(routes,toc,'enumerated',0);
 end
end
function a=pack(routes,elapsed,status,gap)
 a.routes=routes; a.runtime=elapsed; a.status=status; a.gap=gap;
end
function k=measure(I,a,c)
 [C,~,dist]=simulate(I,a.routes,c);
 assert(all(C>0),'Every order must occur once.');
 k.response=C(end)-I.t0; k.tardiness=sum(max(0,C(1:I.n)-I.due));
 k.late=sum(C(1:I.n)>I.due+1e-8); k.distance=dist;
 k.added=max(0,dist-I.baselineDistance); owner=zeros(I.n,1);
 for p=1:I.P, owner(a.routes{p}(a.routes{p}<=I.n))=p; end
 k.reassigned=sum(owner~=I.owner);
 k.deadlineMiss=double(k.response>I.rushDue-I.t0+1e-8);
 k.inversions=0;
 for p=1:I.P
  original=I.base{p}; current=a.routes{p};
   current=current(current<=I.n);
  current=current(I.owner(current)==p);
  for h=1:numel(original)
   ah=find(current==original(h),1);
   if isempty(ah),continue;end
   for j=h+1:numel(original)
    bj=find(current==original(j),1);
    if ~isempty(bj) && ah>bj,k.inversions=k.inversions+1;end
   end
  end
 end
 k.objective=sum(c.weights.*([k.response k.tardiness k.added k.reassigned]./c.reference));
end
function a=optimizeMILP(I,c,incumbent)

 tic; N=I.n+1; P=I.P; arcP=[]; arcI=[]; arcJ=[]; arcd=[];
 for p=1:P
  for j=1:N
   arcP(end+1)=p; arcI(end+1)=0; arcJ(end+1)=j; arcd(end+1)=I.ds(p,j); 
  end
  for i=1:N, for j=1:N
   if i~=j
    arcP(end+1)=p; arcI(end+1)=i; arcJ(end+1)=j;
    arcd(end+1)=I.d(i,j);
   end
  end, end
 end
 E=numel(arcP); ixC=E+(1:N); ixL=E+N+(1:I.n);
 ixZ=E+N+I.n+(1:I.n); ixD=E+N+2*I.n+1; V=ixD;
 obj=zeros(V,1); obj(ixC(end))=c.weights(1)/c.reference(1);
 obj(ixL)=c.weights(2)/c.reference(2);
 obj(ixZ)=c.weights(4)/c.reference(4);
 obj(ixD)=c.weights(3)/c.reference(3);
 lb=zeros(V,1); ub=inf(V,1); ub(1:E)=1; ub(ixZ)=1;

 M=sum(I.service)+N*max([I.d(:);I.ds(:)])/c.speed+1;
 ub(ixC)=M;
 A=sparse(0,V); b=[]; Aeq=sparse(0,V); beq=[];
 for j=1:N 
  v=zeros(1,V); v(find(arcJ==j))=1; Aeq(end+1,:)=sparse(v); beq(end+1)=1;
 end
 for p=1:P
  v=zeros(1,V); v(find(arcP==p & arcI==0))=1;
  A(end+1,:)=sparse(v); b(end+1)=1;
  for i=1:N 
   v=zeros(1,V); v(find(arcP==p & arcI==i))=1;
   v(find(arcP==p & arcJ==i))=-1;
   A(end+1,:)=sparse(v); b(end+1)=0;
  end
 end
 for e=1:E 
  j=arcJ(e); i=arcI(e); v=zeros(1,V);
  v(e)=M; v(ixC(j))=-1;
  if i>0, v(ixC(i))=1; end
  A(end+1,:)=sparse(v);
  b(end+1)=M-I.service(j)-arcd(e)/c.speed;
 end
 for j=1:I.n
  v=zeros(1,V);v(ixC(j))=1;v(ixL(j))=-1;
  A(end+1,:)=sparse(v);b(end+1)=I.due(j);
    v=zeros(1,V);v(find(arcP==I.owner(j) & arcJ==j))=-1;v(ixZ(j))=-1;
  A(end+1,:)=sparse(v);b(end+1)=-1;
 end
 v=zeros(1,V);v(ixZ)=1;A(end+1,:)=sparse(v);b(end+1)=c.maxReassigned;
 v=zeros(1,V);v(ixC(end))=1;A(end+1,:)=sparse(v);b(end+1)=I.rushDue;
 v=zeros(1,V);v(1:E)=arcd;v(ixD)=-1;
 A(end+1,:)=sparse(v);b(end+1)=I.baselineDistance;
 opts=optimoptions('intlinprog','Display','off', ...
  'MaxTime',c.maxSolverSeconds,'RelativeGapTolerance',c.relativeGap);
 try
  [x,~,flag,output]=intlinprog(obj,[1:E ixZ],A,b,Aeq,beq,lb,ub,opts);
 catch ME
  warning('MILP solver error: %s',ME.message); x=[];flag=-99;output=struct();
 end
 if isempty(x)
  a=incumbent;a.runtime=toc;a.status='no incumbent; S1 fallback';a.gap=NaN;return
 end
 routes=cell(P,1);
 for p=1:P
  active=find(x(1:E)>.5 & arcP(:)==p); first=active(arcI(active)==0);
  if isempty(first), continue; end
  j=arcJ(first(1)); routes{p}=j;
  while true
   follow=active(arcI(active)==j); if isempty(follow), break; end
   j=arcJ(follow(1));routes{p}(end+1)=j;
   assert(numel(routes{p})<=N,'Invalid MILP cycle.');
  end
 end
 assert(isequal(sort([routes{:}]),1:N),'Invalid assignment from MILP.');
 if flag==1, status='optimal within gap'; else, status='feasible incumbent'; end
 gap=NaN;if isfield(output,'relativegap'),gap=output.relativegap;end
 a=pack(routes,toc,status,gap);
 if ~strcmp(incumbent.status,'no feasible insertion')
  if measure(I,a,c).objective > measure(I,incumbent,c).objective+1e-8
   a=incumbent; a.runtime=toc; a.status='S1 retained; MILP incumbent worse';
   a.gap=gap;
  end
 end
end
