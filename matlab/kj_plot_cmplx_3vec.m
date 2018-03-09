function [] = kj_plot_cmplx_3vec(E0,E1,E2)

M = numel(E0);

n = M/3;

E0_r = E0(0*n+1:1*n);
E0_t = E0(1*n+1:2*n);
E0_z = E0(2*n+1:3*n);

if exist('E1')
    E1_r = E1(0*n+1:1*n);
    E1_t = E1(1*n+1:2*n);
    E1_z = E1(2*n+1:3*n);
end

if exist('E2')
    E2_r = E2(0*n+1:1*n);
    E2_t = E2(1*n+1:2*n);
    E2_z = E2(2*n+1:3*n);
end

figure();

ax1 = subplot(3,2,1);
plot(ax1,real(E0_r),'Color','k')
hold on
if exist('E1')
    plot(ax1,real(E1_r),'Color','b')
end
if exist('E2')
    plot(ax1,real(E2_r),'Color','r')
end

ax2 = subplot(3,2,2);
plot(ax2,imag(E0_r),'Color','k')
hold on
if exist('E1')
    plot(ax2,imag(E1_r),'Color','b')
end
if exist('E2')
    plot(ax2,imag(E2_r),'Color','r')
end

ax3 = subplot(3,2,3);
plot(ax3,real(E0_t),'Color','k')
hold on
if exist('E1')
    plot(ax3,real(E1_t),'Color','b')
end
if exist('E2')
    plot(ax3,real(E2_t),'Color','r')
end

ax4 = subplot(3,2,4);
plot(ax4,imag(E0_t),'Color','k')
hold on
if exist('E1')
    plot(ax4,imag(E1_t),'Color','b')
end
if exist('E2')
    plot(ax4,imag(E2_t),'Color','r')
end

ax5 = subplot(3,2,5);
plot(ax5,real(E0_z),'Color','k')
hold on
if exist('E1') 
    plot(ax5,real(E1_z),'Color','b')
end
if exist('E2') 
    plot(ax5,real(E2_z),'Color','r')
end

ax6 = subplot(3,2,6);
plot(ax6,imag(E0_z),'Color','k')
hold on
if exist('E1')
    plot(ax6,imag(E1_z),'Color','b')
end
if exist('E2')
    plot(ax6,imag(E2_z),'Color','r')
end

end