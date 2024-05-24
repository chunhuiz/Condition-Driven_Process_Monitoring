%% estimate control limit via KDE
%inputs£º
%       monitoring_index(samples*1)£ºthe control limit calculated using
%               training data.
%       confidence_level£º the confidence level of the control limit, e.g.
%       0.95, 0.99.
%outputs£º
%       Ctrl_limit: control limit for the monitoring index.
function Ctrl_limit=ctrl_limit_compute(monitoring_index,confidence_level)

[f,x]=ksdensity(monitoring_index,'Function','cdf','npoints',1000);
id = find(f-confidence_level>=0);
id = id(1);
Ctrl_limit = x(id);
end