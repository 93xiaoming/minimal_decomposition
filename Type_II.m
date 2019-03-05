%MAIN FUNCTION FOR TYPE II
%Input: parameters of target transformation(theta,psi,phi), and angle
%between two axes (Theta)
%Output: first  arrow: theta of each elementary rotations
%        second arrow: phi of each elementary rotations 

function output = Type_II( theta,psi,phi,Theta )
%check input validity
    if theta<0 ||theta>=pi
        error('invalid theta value');
    end
    if psi<0 ||psi>=pi
        error('invalid psi value');
    end
    if phi<0 ||phi>=4*pi
        error('invalid phi value');
    end
    if Theta<=0 ||Theta>pi/2
        error('invalid Theta value');
    end
    
    [p_even,num_even,beta1] = Bol_even(theta,psi,phi,Theta);
    [p_odd,num_odd]  = Bol_odd(theta,psi,phi,Theta);
    
    pmin = min(p_even, p_odd);   %determine pmin
      
    
    if pmin==1  %pmin==1?
        output = [theta;phi];
        
    elseif pmin==2
        output = p2( theta,psi,phi,Theta );
    elseif pmin==p_odd     %pmin is odd?
        
        output = decomp_odd(theta,psi,phi,Theta,pmin,num_odd);
    elseif pmin==p_even    %pmin is even?
        
        if num_even==1
            ang  = ag( Rot(Theta,0,-beta1)*Rot(theta,psi,phi) );
            opt  = [Theta;beta1];
            opt2 = decomp_odd(ang(1),ang(2),ang(3),Theta,pmin-1,1);
            output = [opt,opt2];
        else
            ang = ag( Rot(theta,psi,phi)*Rot(Theta,0,-beta1) );
            opt  = [Theta;beta1];
            opt2 = decomp_odd(ang(1),ang(2),ang(3),Theta,pmin-1,1);
            output= [opt2,opt];
        end
    end
    output =output +4*pi*(output<0);
end

%SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%minimum even value of p 
function [pmin,num,beta1] = Bol_even( theta,psi,phi,Theta )
    
    A_1 = ( cos(psi)*cos(Theta)*sin(theta)*sin(phi/2)-sin(Theta)*cos(theta)*sin(phi/2) )^2 ;
    A_2 = (sin(psi)*cos(Theta)*sin(theta)*sin(phi/2)-sin(Theta)*cos(phi/2) )^2 ;
    A   = A_1+A_2;
    B   = ( sin(theta)*sin(phi/2) )^2;
    C   = sin(Theta)*sin(theta)*sin(phi/2)* (sin(psi)*sin(phi/2)*cos(theta) - cos(psi)*cos(phi/2)); 
    Lambda =  asin( real( sqrt( (A+B)/2 - sqrt( C^2+ (B-A)^2/4 ) ) ));
    
    A_1_ = ( cos(psi)*cos(Theta)*sin(theta)*sin(-phi/2)-sin(Theta)*cos(theta)*sin(-phi/2) )^2 ;
    A_2_ = (sin(psi)*cos(Theta)*sin(theta)*sin(-phi/2)-sin(Theta)*cos(-phi/2) )^2 ;
    A_   = A_1_+A_2_;
    B_   = ( sin(theta)*sin(-phi/2) )^2;
    C_   = sin(Theta)*sin(theta)*sin(-phi/2)* (sin(psi)*sin(-phi/2)*cos(theta) - cos(psi)*cos(-phi/2));
    Lambda_ = asin( real( sqrt( (A_+B_)/2 - sqrt( C_^2+ (B_-A_)^2/4 ) ) ));
    
    [L,num] = min([Lambda,Lambda_]);
    
    pmin =  2*( ceil(L/Theta)+1 ) ;
    
    if num ==1
        beta1 = pi + angle((B-A)/2+1i*C);
    else
        AA=A_;BB=B_;CC=C_;
        beta1 = -pi - angle((B_-A_)/2+1i*C_);
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%minimum odd value of p 
function [pmin, num, delta] = Bol_odd( theta,psi,phi,Theta)
    delta1  = asin( sin(theta)*sin(phi/2) );
    delta2 = asin( sin(phi/2)*sqrt( (cos(Theta)*cos(psi)*sin(theta)-cos(theta)*sin(Theta))^2+(sin(theta)*sin(psi))^2) );
    
    
    [d,num]= min([abs(delta1), abs(delta2)]);
    
    pmin  = 2*ceil(d/Theta)+1;
    if num==1
        delta=delta1;
    else
        delta=delta2;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%odd piece decomposition of Type II
function output = decomp_odd( theta,psi,phi,Theta,pmin,num )
    
    if num==1
        theta_ = theta;psi_=psi;phi_=phi;
        delta = asin( sin(theta_)*sin(phi_/2) );
        
        gamma   = 2 * asin(sin(delta*2/(pmin-1))/sin(Theta) );
        lambda1 = angle( cos(phi_/2)+1i*sin(phi_/2)*cos(theta_) );
        lambda2 = angle( cos(gamma/2)+1i*sin(gamma/2)*cos(Theta) );
        
        mid=repmat([-2*lambda2,gamma],1,(pmin-3)/2 );
        
        output  = [ repmat([0,Theta],1,(pmin-1)/2 ),0;
                  lambda1-lambda2+psi_, gamma, mid, lambda1-lambda2-psi_];
    else
        ang   = ag(Rot(Theta/2,0,-pi)*Rot(theta,psi,phi)*Rot(Theta/2,0,pi));
        theta_=ang(1);psi_=ang(2);phi_=ang(3);
        delta = asin( sin(theta_)*sin(phi_/2) );
        
        gamma   = 2 * asin(sin(delta*2/(pmin-1))/sin(Theta) );
        lambda1 = angle( cos(phi_/2)+ 1i*sin(phi_/2)*cos(theta_) );
        lambda2 = angle( cos(gamma/2)+1i*sin(gamma/2)*cos(Theta) );
        
        mid=repmat([-2*lambda2,gamma],1,(pmin-3)/2 );
        
        output  = [repmat([Theta,0],1,(pmin-1)/2 ),Theta;
                  lambda1-lambda2+psi_, gamma, mid, lambda1-lambda2-psi_];
        
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input: transformation matrix 
%Output: corresponding parameters (theta,psi,phi)

function ag = ag( R )     
%azimuthal angle    
    if imag(R(2,1)*i)>0        
        psi=angle(R(2,1)*i);
    elseif imag(R(2,1)*i)<0
        psi=angle(-R(2,1)*i);
    else 
        psi=0;
    end
%rotation angle
    phi=2*acos(real(R(1,1))); 
    imag(R(2,1)*i);
    if imag(R(2,1)*i)<0
        phi=4*pi-phi;
    end
%polar angle
    if abs(sin(phi/2))<10e-8
        theta=0;
    else 
        theta=acos((-imag(R(1,1))/sin(phi/2)));
    end

    ag=[theta psi phi];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input: target transformation parameters
%Output: transformation matrix
function Rot = Rot(theta,psi,phi)   
sx=[0 1;1 0];sy=[0 -i;i 0];sz=[1 0;0 -1];
Rot=expm(-i*(sin(theta)*cos(psi)*sx+sin(theta)*sin(psi)*sy+cos(theta)*sz)*phi/2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output: decomposition for p=2
function output = p2( theta,psi,phi,Theta )
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
    end
end  







