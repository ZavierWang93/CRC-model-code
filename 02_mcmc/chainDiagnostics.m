function report = chainDiagnostics(draws)

[n,m,p] = size(draws);
if m < 2 || n < 20
    error('MeanageDiagnostics:Chains','Use at least two independent chains with >=20 draws.');
end
half = floor(n/2);
split = cat(2,draws(1:half,:,:),draws(end-half+1:end,:,:));
rhat = nan(1,p); bulk = zeros(1,p); tail = zeros(1,p);
for k = 1:p
    x = split(:,:,k);
    if any(~isfinite(x),'all') || any(all(x==x(1,:),1))
        continue;
    end
    z = rankNormalize(x);
    folded = rankNormalize(abs(x-median(x,'all')));
    rhat(k) = max(splitRhat(z),splitRhat(folded));
    bulk(k) = multiESS(z);
    q = prctile(x(:),[5 95]);
    tail(k) = min(multiESS(double(x<=q(1))),multiESS(double(x<=q(2))));
end
report = struct('rhat',rhat,'ess_bulk',bulk,'ess_tail',tail, ...
    'max_rhat',max(rhat),'min_bulk_ess',min(bulk),'min_tail_ess',min(tail), ...
    'draws_per_chain',n,'num_chains',m, ...
    'passed',all(isfinite(rhat) & rhat<1.01 & bulk>=max(400,100*m) & tail>=max(400,100*m)));
if any(~isfinite(rhat)), report.max_rhat=Inf; end
end

function z = rankNormalize(x)
ranks = tiedrank(x(:));
z = reshape(norminv((ranks-3/8)/(numel(ranks)+1/4)),size(x));
end

function r = splitRhat(x)
n = size(x,1);
W = mean(var(x,0,1));
B = n*var(mean(x,1),0,2);
if W<=0, r=NaN; return; end
r = sqrt(((n-1)/n*W+B/n)/W);
end

function ess = multiESS(x)
[n,m] = size(x);
variances = var(x,0,1);
W = mean(variances);
varPlus = (n-1)/n*W + var(mean(x,1),0,2);
if ~isfinite(varPlus) || varPlus<=0
    ess=0; return;
end
centered = x-mean(x,1);
spectrum = fft(centered,2^nextpow2(2*n-1),1);
acov = ifft(abs(spectrum).^2,[],1,'symmetric');
acov = acov(1:n,:)/n;
rho = 1-(W-mean(acov,2))/varPlus;
rho(1) = 1;
pairs = rho(1:2:end-1)+rho(2:2:end);
stop = find(pairs<=0 | ~isfinite(pairs),1);
if ~isempty(stop), pairs=pairs(1:stop-1); end
if isempty(pairs), ess=0; return; end
pairs = cummin(pairs);
tau = max(1,-1+2*sum(pairs));
ess = min(n*m,n*m/tau);
end
