clc
clear
n = 1000
Tmin = 1
Tmax = 500
T = linspace(Tmin,Tmax,n)
Na = 6.022D23
Kb = 1.3806D-23
θ_D = [400 428 470 630] //[Mg Al Fe Cr]
θ_E = 0.806*θ_D
//Dulong Petit Law
Cv_DP = 3*Na*Kb
show_window(1)
plot(T,Cv_DP,'--r')
//Einstein Theory
for i=1:n
    for j=1:4
    Cv_E(i,j) =3*Na*Kb*((θ_E(j)/T(i))^2)*(exp(θ_E(j)/T(i))/(((exp(θ_E(j)/T(i)))-1)^2))
    end
end
show_window(2)
for k=1:4
    //T_E(k)=T(1:4)/θ_E(k)
    subplot(2,2,k)
    plot((T'/θ_E(k)),Cv_E(:,k))
end
//Debay Model
for ii=1:4
    x = θ_D(ii)./T
end
function y=f(x)
    y=(x^4*exp(x))/(((exp(x))-1)^2)
endfunction
for i = 1:n
    for j=1:4
        Z(i,j)=intg(0,(θ_D(j)/T(i)),f)
        Cv_D(i,j)=9*Na*Kb*((T(i)/θ_D(j))^3)*Z(i,j)
     end
end
show_window(3)
for k=1:4
    //T_E(k)=T(1:4)/θ_E(k)
    subplot(2,2,k)
    plot((T'/θ_D(k)),Cv_D(:,k))
end
//lower temperature in Debay
for i=1:n
    for j=1:4
        Cv_DL(i,j)=((12*(%pi)^4)/5)*Na*Kb*((T(i)/θ_D(j))^3)
        Cv_DPL(i,j) = 3*Na*Kb
    end
end
show_window(4)
for k=1:4
    //T_E(k)=T(1:4)/θ_E(k)
    subplot(2,2,k)
    plot((T^3)',Cv_DL(:,k),'*r')
    plot((T^3)',Cv_D(:,k),'-k','linewidth',1.5)
    at=gca(); 
    at.data_bounds=[0,0;8e6,25]
end
//Comparison of three laws----------------
show_window(5)
for k=1:4
    subplot(2,2,k)
    plot((T^3)',Cv_D(:,k),'-k')
    plot((T^3)',Cv_DL(:,k),'-r')
    plot((T^3)',Cv_E(:,k),'-o')
    at=gca(); 
    at.data_bounds=[0,0;1.3e6,25]
end
show_window(6)
for k=1:4
    subplot(2,2,k)
    plot(T',Cv_E(:,k),'-o')
    plot(T',Cv_D(:,k),'-k','linewidth',1.5)
    plot(T',Cv_DPL(:,k),'--r')
    //at=gca(); 
    //at.data_bounds=[0,0;1.3e6,25]
end
