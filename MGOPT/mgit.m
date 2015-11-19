% Perform one iteration of multigrid.
% Can only be used after using mg.m or fmg.m
%----------------------------------------------------------------------
global GRAPH_N_OLD GRAPH_INDEX_OLD it W_level WS
%----------------------------------------------------------------------
% Update the multi-grid solution: given current solution v
%----------------------------------------------------------------------
h_MG_graph = findobj('Tag', 'multigrid_graph');
if (~isempty(h_MG_graph))
  GRAPH_INDEX_OLD = GRAPH_INDEX_OLD+8;
  figure(h_MG_graph);
  hold on;
else
  GRAPH_N_OLD = N(1);
  GRAPH_INDEX_OLD = 1;
  h_MG_graph = figure('Tag', 'multigrid_graph', ...
					  'IntegerHandle', 'off', ...
					  'Name', 'Multigrid graph', ...
					  'Numbertitle', 'off', ...
					  'Visible', 'off');
end

title('Progress of Multigrid Iteration');
set(gca, 'XTickLabel', [], ...
		 'YTickMode', 'manual', 'YTick', [1:numel(N)], 'YTickLabel', N, ...
		 'YLimMode', 'manual', 'YLim', [1, numel(N)]);
hold on;
%----------------------------------------------------------------------
more off;
it = it + 1;
step_bnd = 0;
v  = mgrid(v,fnl,0,step_bnd);
% doplot(it,v,W_level);
report_results(N);
MGiter=0;
for i=1:size(NF,2),MGiter=MGiter+sum(NF(2:3,i),1)/(4^(i-1));end
[f,g]=sfun(v);
ResMg(it+1,1:3)=[MGiter,f,norm(v-WS(:))];
more on;
%----------------------------------------------------------------------
figure(findobj('Tag', 'multigrid_graph'));
hold off;
