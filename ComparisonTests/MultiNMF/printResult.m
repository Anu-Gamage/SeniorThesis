function [acc_val, nmi_val] = printResult(X, label, K, kmeansFlag)
if(~exist('kmeansFlag','var'))
		kmeansFlag = 1;
	end
    for i=1:1
        if kmeansFlag == 1
            indic = litekmeans(X, K, 'Replicates',20);
        else
            [~, indic] = max(X, [] ,2);
        end
        result = bestMap(label, indic);
        [ac(i), nmi_value(i), cnt(i)] = CalcMetrics(label, indic);
    end
  %  disp(sprintf('ac: %0.3f\t%d/%d\tnmi:%0.3f\t', mean(ac), round(mean(cnt)), length(label), mean(nmi_value)));
    acc_val = mean(ac);
    nmi_val = mean(nmi_value);
end
