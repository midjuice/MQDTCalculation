close all;
clear;
clc;

gJ = 2.0023;

gI =5e-4;
Ias = [1, 1.5, 2.5, 4];
Ias_print = ["1", "3/2", "5/2", "4"];
Ehf = 300;

left = .1;
right = 200;
points = 5000;
Bs = linspace(left, right, points);
output = zeros(1, points);

t = tiledlayout(2,2,"TileSpacing","none");
for s=1:length(Ias)
    Ia = Ias(s);
    mfs = -Ia-1/2:1:Ia+1/2;
    % subplot(2,2,s);
    nexttile;
    draw_structure(points, output, Bs, mfs, gJ, gI, Ia, Ehf)
    xlim([0,240]);
    set(gca,'xticklabel',[],'yticklabel',[]);
    xL=xlim;
    yL=ylim;
    text(10, 0.95*yL(2), "I="+Ias_print(s), 'HorizontalAlignment','left','VerticalAlignment','top')
end
xlabel(t, 'x');
ylabel(t, 'Zeeman Energy');

% saveas(t, "hyperfine_structure.jpg");
% print(gcf,'hyperfine_structure.png','-dpng','-r300'); 
% https://stackoverflow.com/questions/32789901/how-to-save-a-high-resolution-figure-in-matlab

function draw_structure(points, output, Bs, mfs, gJ, gI, Ia, Ehf)
    for count=1:length(mfs)
        mf = mfs(count);
        if (-Ia-1/2<mf) && (mf<Ia+1/2)
            alpha = 1;
            draw(points, output, Bs, mf, alpha, gJ, gI, Ia, Ehf)
            alpha = -1;
            draw(points, output, Bs, mf, alpha, gJ, gI, Ia, Ehf)
        else
            alpha = 0;
            draw(points, output, Bs, mf, alpha, gJ, gI, Ia, Ehf)
        end
    end
end

function draw(points, output, Bs, mf, alpha, gJ, gI, Ia, Ehf)
    for i=1:points
        B = Bs(i);
        output(i) = Ethfunc(B,mf,alpha,gJ,gI,Ia,Ehf);
    end
    plot(Bs, output);
    if mod(Ia, 1)~=0
        text(Bs(end), output(end), sprintf('$$(%2.0f)_{%d}$$',mf,alpha), 'Interpreter', 'latex', 'FontSize', 7);
    else
        text(Bs(end), output(end), sprintf('$$(%2.1f)_{%d}$$',mf,alpha), 'Interpreter', 'latex', 'FontSize', 7);
    end
    hold on;
end
