maturity = 3; % year
rfr = 0.06;
spot = 1;
strike = 1.1;
dt = 1;

sqrt_dt = sqrt(dt);
time = dt:dt:maturity;

for divYield = [0];
    
    %       Results from regular Monte-Carlo
    S = [1.09 1.16 1.22 0.93 1.11 0.76 0.92 0.88;
        1.08 1.26 1.07 0.97 1.56 0.77 0.84 1.22;
        1.34 1.54 1.03 0.92 1.52 0.9 1.01 1.34]';%?
    %     S = exp(logS);
    ECF = max(strike-S,0); % ECF = Exercise Cash Flow %?
    Moneyness = ECF>0;
    PVCCF = zeros(size(S,1),size(S,2)-1); % PVCCF = Present Value of Continuing Cash Flow
    CondExpect = zeros(size(S,1),size(S,2)-1);
    
    %       Initialization for LSMC
    CFM = zeros(size(S)); % CFM = Cash Flow Matrix
    CFM(:,end) = ECF(:,end);
    
    for maxSpower = 2
        for i=length(time)-1:-1:1
            disp(['Time=',num2str(i)]);
            PVCCF(:,i) = CFM(:,i+1)*exp(-rfr*dt);
            for j=i+2:length(time)
                PVCCF(:,i) = PVCCF(:,i) + CFM(:,j)*exp(-rfr*(j-i)*dt);
            end
            [inMoneyPaths,~,~] = find(ECF(:,i));
            basisVec = ones(length(inMoneyPaths),1);
            for j=1:maxSpower
                basisVec = [basisVec , S(inMoneyPaths,i).^j];
            end
            alphas = basisVec\PVCCF(inMoneyPaths,i)
            CondExpect(inMoneyPaths,i) = basisVec*alphas;
            CFM(:,i) = (CondExpect(:,i)<ECF(:,i)).*ECF(:,i);
            for j=i+1:size(CFM,2)
                CFM(:,j) = (CondExpect(:,i)>=ECF(:,i)).*CFM(:,j);
            end
        end
        BermudanPut = mean(CFM*exp(-rfr*time'));
    end
end

ECF
PVCCF
CondExpect
CFM

disp(['BermudanPut=',num2str(BermudanPut)]);
