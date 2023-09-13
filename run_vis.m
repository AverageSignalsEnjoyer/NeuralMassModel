clear all
close all

load y

f = figure;
f.Position = [1138 689 862 649];
set(f, 'color', 'w')

hold on
box on
zlim([-1 1])
view(30, 30)
% dim = [6 6];
colormap('hsv');
caxis([-0.5 0.5])
zlabel('Membrane depolarization (mV)')
xticks([0:5]*0.4)
xticklabels(string([1:6]))
yticks([0:5]*0.4)
yticklabels(string([1:6]))
xlabel('Neural Mass Index (x-dim)')
ylabel('Neural Mass Index (y-dim)')
ax = gca;
ax.FontSize = 16;

vidfile = VideoWriter('without_feedback.mp4','MPEG-4');
vidFile.Quality = 100;
open(vidfile);

[X,Y,Z] = sphere(50);
X = X*0.2;
Y = Y*0.2;
Z = Z*0.2;

% ---- Init grid ----
objs = {};
for i = 1:dim(1)
    for j = 1:dim(2)

        idx = sub2ind(dim, i, j);
        objs{idx} = surf(X+(i-1)*0.4,Y+(j-1)*0.4,Z);
        objs{idx}.EdgeColor = 'none';

    end
end

for t = 1:size(y, 2)

    for i = 1:dim(1)
        for j = 1:dim(2)

            idx = sub2ind(dim, i, j);
            objs{idx}.ZData = Z+y(idx, t)/12;

        end
    end

    drawnow
    F(t) = getframe(gcf);
    writeVideo(vidfile,F(t));
end
close(vidfile)

