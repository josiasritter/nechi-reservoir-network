function [fdc_deficit,fdc_surplus] = FDC_Shift(dCurve_Pre,dCurve_Post,pEmp,percentileRanges)

    [r] = max(size(percentileRanges));

    fdc_deficit = zeros(r-1,1);
    fdc_surplus = zeros(r-1,1);

    diff = dCurve_Post - dCurve_Pre;
    deficits = diff .* (diff < 0);
    surplus = diff .* (diff > 0);
   
    for i = 1:r-1
    
        a1 = pEmp >= percentileRanges(i);
        a2 = pEmp <= percentileRanges(i+1);
        
        a0 = sum(a1 .* a2 .* dCurve_Pre);
        
        fdc_deficit(i) = -1 * sum(a1 .* a2 .* deficits) / a0; 
        fdc_surplus(i) = sum(a1 .* a2 .* surplus)/ a0; 
        
    end
    
end

