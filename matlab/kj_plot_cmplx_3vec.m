function [] = kj_plot_cmplx_3vec(E0,E1,E2)

M = numel(E0);

n = M/3;

E0_1 = E0(0*n+1:1*n);
E0_2 = E0(1*n+1:2*n);
E0_3 = E0(2*n+1:3*n);

if exist('E1')
    E1_1 = E1(0*n+1:1*n);
    E1_2 = E1(1*n+1:2*n);
    E1_3 = E1(2*n+1:3*n);
end

if exist('E2')
    E2_1 = E2(0*n+1:1*n);
    E2_2 = E2(1*n+1:2*n);
    E2_3 = E2(2*n+1:3*n);
end

figure();

ax1 = subplot(3,2,1);
plot(ax1,real(E0_1),'Color','k')
hold on
if exist('E1')
    plot(ax1,real(E1_1),'Color','b')
end
if exist('E2')
    plot(ax1,real(E2_1),'Color','r')
end

ax2 = subplot(3,2,2);
plot(ax2,imag(E0_1),'Color','k')
hold on
if exist('E1')
    plot(ax2,imag(E1_1),'Color','b')
end
if exist('E2')
    plot(ax2,imag(E2_1),'Color','r')
end

ax3 = subplot(3,2,3);
plot(ax3,real(E0_2),'Color','k')
hold on
if exist('E1')
    plot(ax3,real(E1_2),'Color','b')
end
if exist('E2')
    plot(ax3,real(E2_2),'Color','r')
end

ax4 = subplot(3,2,4);
plot(ax4,imag(E0_2),'Color','k')
hold on
if exist('E1')
    plot(ax4,imag(E1_2),'Color','b')
end
if exist('E2')
    plot(ax4,imag(E2_2),'Color','r')
end

ax5 = subplot(3,2,5);
plot(ax5,real(E0_3),'Color','k')
hold on
if exist('E1') 
    plot(ax5,real(E1_3),'Color','b')
end
if exist('E2') 
    plot(ax5,real(E2_3),'Color','r')
end

ax6 = subplot(3,2,6);
plot(ax6,imag(E0_3),'Color','k')
hold on
if exist('E1')
    plot(ax6,imag(E1_3),'Color','b')
end
if exist('E2')
    plot(ax6,imag(E2_3),'Color','r')
end

end