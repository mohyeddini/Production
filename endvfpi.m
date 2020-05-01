%%%%%%%%%%%%%%%%%%%%%%%%%%%%AmirMohyeddini%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d=4.5/12;%ft
WOR=0;
GOR=1630;%scf/stb
GLR=1630;%scf/bbl
mug=0.02;%cp
muw=1;%cp
API=39;
SGw=1;
SGg=0.75;
gamag=0.05;
pbar=14.7;
qo=1380;%bbl/day
qw=0;%bbl/day
ql=qo+qw;
f=0.006;
sigmaw=0;%dyne/cm
sigmao=40;%dyne/cm
%f=input('please Enter f From moody for NRe and epsilon/d:');
%%%%%%%%%%%%%%%
deltapo=10;%psi
pold=5000;
epsilon=10;
deltahn=[];
P=[];
depth=[];
deltahf=zeros(1,200);
depthh=zeros(1,200);
for ff=1:200
   deltahf(ff)=70; 
end
for ff=13930:-70:0
   depthh=ff;
   depth=[depth,depthh];%#ok
end

for hh=1:1:200
    deltah=deltahf(hh);
    deltapo=10;
    epsilon=10;
    while epsilon>=0.001

    T=temp(depth(hh));
    pnew=pold-deltapo;

    % Specific Gravity of the oil
    SGo=141.5/(131.5+API);

    % Find the weight associated with 1 bbl of stock-tank liquid
    %m=SGo*(350)*(1/(1+WOR))+SGo*(350)*(WOR/(1+WOR))+(0.0764)*(GLR)*(SGg);

    % specific weight of the liquid phase
    gamal=62.4*(SGo*(1/(1+WOR))+SGw*(WOR/(1+WOR)));%lb/ft^3
    gamao=62.4*SGo;

    % Zave for gas
    %Zave=f(Tr,Pr);
    Zave=0.85;

    % Average specific weight of the gas phase
    % gamag=SGg*(0.0764)*(Pave/14.7)*(520/5)*(1/Zave);

    %Find Rs at Tave and Pave
    Rs=SGg*((pnew/18)*(10^(0.0125*(API)))/(10^(0.00091*(T+460))))^1.2048;
    %Rs=10^(0.3818-5.506*log10(So)+2.902*log10(Sgs)+1.327*log10(Ps)-0.7355*log10(Ts));
    %So is the stock tank oil specific gravity,Sgs, Ps and Ts, are
    %the separator gas specific gravity, pressure, and temperature


    % %rhoo
    % rhoal=38.52*10^(-0.00326*API)+(94.75-33.93*log10(API))*log10(Sg);
    % rhoao=((Rs/380)*(28.96*Sg)+5.61*62.4*So)/(5.61+((Rs/380)*(28.96*Sg)/rhoal));
    % deltarhop(i)=(0.167+16.181*10^(-0.0425*rhoao))*(p(i)/1000)-0.01*(0.299+263*10^(-0.0603*rhoao))...
    %     *((p(i)/1000)^2);
    % deltarhot(i)=(0.0133+152.4*(rhoao+deltarhop(i))^-2.45)*(T(i)-60)-(8.1*(10^-6)-0.0622*...
    %     10^-0.0764*(rhoao+deltarhop(i)))*(T(i)-60)^2;
    % 
    % rhoo(i)=rhoao+deltarhop(i)-deltarhot(i);


    %average viscosity of the oil from correlation
    x=(T^-1.163)*exp(6.9824-0.04658*API);
    muod=(10^x)-1;
    AA=10.715*(Rs+100)^-0.515;
    B=5.44*(Rs+150)^-0.338;
    muo=AA*(muod^B);

    % %muo Beggs and Robinson
    % AA=10^(3.0324-0.02023*API-1.163*log10(T(i)));
    % muod=(10^AA)-1;
    % CC=10.715*(Rs+100)^-0.515;
    % BB=5.44*(Rs+150)^-0.338;
    % muob=(CC*muod)^BB;
    % DD(i)=2.6*(p(i)^1.187)*exp(-11.513-8.98*(10^-5)*p(i));
    % muo(i)=muob*(p(i)/pb)^DD(i);

    % % gas
    % % mug
    % aaa(i)=((9.379+0.0160*M)*T(i)^1.5)/(209.2+19.26*M+T(i));
    % bbb(i)=3.448+0.01009*M+(986.4/T(i));
    % ccc(i)=2.4-0.2*bbb(i);
    % mug(i)=(10^-4)*aaa(i)*exp(bbb(i)*(rhog/62.43)^ccc(i));

    %Liquid mixture viscosity
    mul=muo*(1/(1+WOR))+muw*(WOR/(1+WOR));
    % mu1=muh+landan2+landaco2+landah2s;
    % landan2=(yn2*(10^-3))*(9.59+8.48*log10(Sg));
    % landaco2=(yco2*(10^-3))*(6.24+9.08*log10(Sg));
    % landah2s=(yh2s*(10^-3))*(3.73+8.49*log10(sg));

    %Liquid mixure surface tension
    sigmal=sigmao*(1/(1+WOR))+sigmaw*(WOR/(1+WOR));
    %**************************************************************************
    % Find Bo and Tave
    % F=Rs*sqrt(gamag/gamao)+(1.25/100)*(Tave);
    % Bo=0.972+0.00147*F^1.175;
    %BO Standig correlation
    Bo=0.9759+0.000120*(Rs*sqrt(SGg/SGo)+1.25*T)^1.2;

    % %Pb Standing
    % aa=0.00091*T(i)-0.0125*(API);
    % pb=18.2*(((Rs/Sg)^0.83)*(10^aa)-1.4);

    % %Co
    % Co(i)=(-1433+5*Rs+17.2*T(i)-1180*Sg+12.61*API)/((10^5)*p(i));

    %Find the turbin flow area
    Ap=(3.14*(d^2))/4;

    Vsl=((5.615*ql)/(86400*Ap))*(Bo*(1/(1+WOR))+(WOR/(1+WOR)));%ft/s
    Vsg=((ql*(GLR-Rs*(1/(1+WOR))))/(86400*Ap))*(14.7/pnew)*((T+460)/520)*(Zave/1);%ft/s
    rhol=gamal;
    rhog=gamag;
    g=32.2;
    gc=32.2;
    Vm=Vsl+Vsg;
    %%%%%%%%%%%%%%%%%%%%%%%%%%MAIN CALCULATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Nx=Vsg*((rhog/0.0764)^(1/3))*((72*rhol)/(62.4*rhol))^(1/4);
    Ny=Vsl*((72*rhol/62.4*sigmal))^(1/4);

    N1=0.51*(100*Ny)^0.172;
    N2=8.6+3.8*Ny;
    N3=70*(100*Ny)^(-0.152);

    %Vsl=superficial liquid velocity,ft/sec
    %Vsg=superficial gas velocity,ft/sec
    %rhog=gas density .lbm/cu ft
    %rhol=loquid density,lbm/cu ft
    %sigmal=gas-liquid interfacial tension ,dynes/cm

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%Bubble Flow%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if Nx<N1

        Vbs=(1.41*((sigmal*g*(rhol-rhog))/rhol^2)^(1/4));
        Vbf=1.2*Vm+Vbs;%Vbf=bubble rise velocity in the flowing stream
        Hl=1-(Vsg/Vbf);
        NRe=(1488*rhol*Vm*d)/mul;
        disp(NRe);
        dppdzf=(f*rhol*Hl*(Vm^2))/(2*gc*d/12);
        dppdze=(g/gc)*(rhol*Hl+rhog*(1-Hl));
        dppdz=dppdzf+dppdze;
        disp('Bubble');

    %correct%%%%%%%%%%%%%%%%%%%%%%%%%%%Slug Flow%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    elseif ((N1<Nx<N2)&&Ny<4)||((N1<Nx<26.5)&&Ny>=4)%#ok  

        Nv=real(1488*(sqrt((d^3)*g*rhol*(rhol-rhog)))/mul);

        %calculate m
        if Nv>=250
            m=10;
        elseif 250>Nv>18 %#ok
            m=69*Nv^-0.35;
        elseif Nv<=18
            m=25;
        end


        NE=(14610*(d^2)*(rhol-rhog))/sigmal;
        c=real(0.345*(1-exp(-0.029*Nv))*(1-exp((3.37-NE)/m)));

        Vbs=real(c*sqrt((g*d*(rhol-rhog)/rhol)));
        Vbf=1.2*Vm+Vbs;
        Hl=1-(Vsg/Vbf);
        NRe=(1488*rhol*Vm*d)/mul;
        disp(NRe);
        dppdze=(g/gc)*(rhol*Hl+rhog*(1-Hl));
        dppdzf=(f*rhol*Hl*(Vm^2))/(2*gc*d);
        dppdz=dppdzf+dppdze;
        disp('slug');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Mist Flow%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    elseif (Nx>N3&&Ny<4)||(Nx>26.5&&Ny>4)

          rhon=rhol*(Vsl/Vm)+rhog*(Vsg/Vm);
    %     NRe=(rhog*Vsg*d)/mug;
    %     Nwe=(rhog*(Vsg^2)*epsilon)/sigmal;%weber number
    %     Nmu=(mul^2)/(rhol*sigmal*epsilon);
    %  
    %     if (Nwe*Nmu)<0.005
    %          epsilonpd=(0.0749*sigmal)/(rhog*Vsg^2*d);
    %     elseif (Nwe*Nmu)>0.005
    %          epsilonpd=((0.3713*sigmal)/(rhog*(Vsg^2)*d))*(Nwe*Nmu)^0.302;
    %     end
    %     
    %     if epsilonpd>0.05
    %          f=((1/(4*log10(0.27*epsilonpd))^2)+0.067*(epsilonpd)^1.73)*4;
    %          disp(f);
    %     else
    %         disp(NRe);
    %         f=input('please Enter f From moody for NRe and epsilon/d:');    
    %     end
    %     
          NRe=(1488*rhol*Vsl*d)/mul;
          disp(NRe);
          dppdzf=(f*gc*(vsg)^2)/(2*gc*d);
    %     s=input('please inter s for Mist Flow');
    %     Vs=s/((rhol/(sigmal*g))^(1/4));
    %     Hl=(Vs-Vm+sqrt(((Vm-Vs)^2)+4*Vs*Vsl))/(2*Vs);

        s=0;
        Hl=(1/(1+Vsg/Vsl));
        rhos=rhol*Hl+rhog*(1-Hl);
        dppdze=(g/gc)*rhos;

        %acceleration Term
        dppdza=((Vm*Vsg*rhon)/(gc*Pnew))*dppdz;
        Ek=(Vm*Vsg*rhon)/(gc*Pnew);

        dppdz=((dppdzel)+(dppdzf))/(1-Ek);
        disp('Mist Flow');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%Transient Flow%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif (N2<Nx<N3)&&Ny<4 %#ok


        %%%%%%%%%%%slug%%%%%%%%%%%%
        Nv=real(1488*(sqrt((d^3)*g*rhol*(rhol-rhog)))/mul);

        %calculate m
        if Nv>=250
            m=10;
        elseif 250>Nv>18 %#ok
            m=69*Nv^-0.35;
        elseif Nv<=18
            m=25;
        end


        NE=(14610*(d^2)*(rhol-rhog))/sigmal;
        c=real(0.345*(1-exp(-0.029*Nv))*(1-exp((3.37-NE)/m)));

        Vbs=real(c*sqrt((g*d*(rhol-rhog)/rhol)));
        Vbf=1.2*Vm+Vbs;
        Hl=1-(Vsg/Vbf);
        NRe=(1488*rhol*Vm*d)/mul;
        disp(NRe);
        dppdze=(g/gc)*(rhol*Hl+rhog*(1-Hl));
        dppdzf=(f*rhol*Hl*(Vm^2))/(2*gc*d);
        dppdzs=dppdzf+dppdze;
        disp('slug');


        %%%%%%%%%%%Mist%%%%%%%%%%%%
         rhon=rhol*(Vsl/Vm)+rhog*(Vsg/Vm);
    %     NRe=(rhog*Vsg*d)/mug;
    %     Nwe=(rhog*(Vsg^2)*epsilon)/sigmal;%weber number
    %     Nmu=(mul^2)/(rhol*sigmal*epsilon);
    %  
    %     if (Nwe*Nmu)<0.005
    %          epsilonpd=(0.0749*sigmal)/(rhog*Vsg^2*d);
    %     elseif (Nwe*Nmu)>0.005
    %          epsilonpd=((0.3713*sigmal)/(rhog*(Vsg^2)*d))*(Nwe*Nmu)^0.302;
    %     end
    %     
    %     if epsilonpd>0.05
    %          f=((1/(4*log10(0.27*epsilonpd))^2)+0.067*(epsilonpd)^1.73)*4;
    %          disp(f);
    %     else
    %         disp(NRe);
    %         f=input('please Enter f From moody for NRe and epsilon/d:');    
    %     end
    %     
          NRe=(1488*rhol*Vsl*d)/mul;
          disp(NRe);
          dppdzf=(f*gc*(vsg)^2)/(2*gc*d);
    %     s=input('please inter s for Mist Flow');
    %     Vs=s/((rhol/(sigmal*g))^(1/4));
    %     Hl=(Vs-Vm+sqrt(((Vm-Vs)^2)+4*Vs*Vsl))/(2*Vs);

        s=0;
        Hl=(1/(1+Vsg/Vsl));
        rhos=rhol*Hl+rhog*(1-Hl);
        dppdze=(g/gc)*rhos;

        %acceleration Term
        dppdza=((Vm*Vsg*rhon)/(gc*Pnew))*dppdz;
        Ek=(Vm*Vsg*rhon)/(gc*Pnew);

        dppdzm=((dppdzel)+(dppdzf))/(1-Ek);
        disp('Mist Flow');

     %%%%%%%%%%%Transient%%%%%%%%%%%%%

        A=(N3-Nx)/(N3-N2);
        B=(Nx-N2)/(N3-N2);
        dppdz=A*(dppdzs)+B*(dppdzm);
        disp(dppdz);
        disp('Transient');   



    end
    dppdz=dppdz/144;
    deltapn=dppdz*deltah;

    epsilon=abs(deltapn-deltapo);
    deltapo=deltapn;
    pnew=pold-deltapn;
    

    end
    pold=pnew;
    P=[P,pnew];%#ok
end

%Result
disp(P);
plot(P,-depth);

%Temperature
function y=temp(deltah)
y=-3E-07*deltah^2+0.0101*deltah+56.196;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%