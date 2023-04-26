function [EXIT_mach] = EXIT_mach_calc(AeAt,gam0)
    % Generate mach_table
    resolution=500;
    mach_table=linspace(1,6,resolution);
    
    for i=1:length(mach_table)
        AeAt_table(i)=1/mach_table(i)*((1+0.5*(gam0-1)*mach_table(i)^2)/(1+0.5*(gam0-1)))^((gam0+1)/(2*(gam0-1)));
        if(AeAt_table(i)<=AeAt)
            EXIT_mach=mach_table(i);
        end
    end
end