%MAIN FUNCTION FOR TYPE I
%Input: parameters of target transformation(theta,psi,phi), and axes regime
%(defined by Theta)
%Output: first  arrow: theta of each elementary rotations
%        second arrow: phi of each elementary rotations 

function output = Type_I( theta,psi,phi,Theta )
% check input validity
    if theta<0 ||theta>=pi
        error('invalid theta value');
    end
    if psi<0 ||psi>=pi
        error('invalid psi value');
    end
    if phi<0 ||phi>=4*pi
        error('invalid phi value');
    end
    if Theta<=0 ||Theta>pi
        error('invalid Theta value');
    end

    p=1;  %pmin==1?
    if theta<Theta && psi==0
        output=[theta;phi];
    elseif phi==0 || phi==2*pi
        output=[theta;phi];
    elseif theta==0;
        output=[0;phi];
    else
        p=2;
    end
    
    if p==2  %pmin==2?
        [bol,opt2] = p2(theta,psi,phi,Theta);
        if bol==1
            output = opt2;
        elseif Theta>=pi/2  %pmin==3?
            output=Type_II(theta,psi,phi,pi/2);
        else                %pmin>=3
            output=Type_II(theta,psi,phi,Theta);
        end
      
    
    end
    output =output +4*pi*(output<0);  
    
%SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%p=2 decomposition    
function [ bol,output ] = p2( theta,psi,phi,Theta )
    a=( sin(psi)*cos(phi/2)+cos(psi)*sin(phi/2)*cos(theta))/(sin(phi/2)*sin(theta));
    b=(-sin(psi)*cos(phi/2)+cos(psi)*sin(phi/2)*cos(theta))/(sin(phi/2)*sin(theta));

  
    if a>=cot(Theta);
        bol=1;
        
        phi1=2*pi+ (2*acos( cos(phi/2)*cos(psi)-sin(phi/2)*sin(psi)*cos(theta))-2*pi)*sign(sin(phi/2));
        phi2=mod(-2*psi,4*pi);
        
        output = [acot(a)+pi*(acot(a)<0), 0;phi1, phi2];
    elseif b>=cot(Theta);
        bol=1;
        
        phi1=2*psi;
        phi2=2*pi+ (2*acos( cos(phi/2)*cos(psi)+sin(phi/2)*sin(psi)*cos(theta))-2*pi)*sign(sin(phi/2));
                
        
        output = [0, acot(b)+pi*(acot(b)<0);phi1, phi2];
    else 
        bol=0;
        output=[0;0];

end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        
        
        