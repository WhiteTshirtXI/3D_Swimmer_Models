clear
figure
m = 1;
A = 1:30;
gridSize = 256;
for t= 0.1:0.1:3
    filename = sprintf('./data/imworm_3D_VE_n%d_t%f.mat',gridSize,t);
    load(filename);
    Sh = Shat;
    for i = 1:6
        S1(:,:,:,i) = real(ifftn(Sh(:,:,:,i)));
    end
    n = sqrt(sum(S1(:).^2));
    A(m) = n;
    plot(t,n,'bo');
    text(t,n,num2str(n));
    m = m+1;
    hold on
end
tm = 0.1:0.1:3;
line(tm, A);
titl = sprintf('Stress Norm on %d grid', gridSize);
title(titl);
xlabel('Time - Seconds');
ylabel('L2 Norm');
hold off
