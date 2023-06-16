function [MCW]= MultiChannelWiener(x,a,SignalPower,Rn)

    
    MVDR = MVDR(x,a) 
    MCW = ((SignalPower^2)*MVDR)/(SignalPower^2 +inv(a'*inv(Rn)*a))


    

