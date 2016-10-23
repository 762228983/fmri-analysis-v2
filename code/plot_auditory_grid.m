function G_plot = plot_auditory_grid(G, varargin)

X = [G.grid_data{1}(:); G.grid_data{2}(:)];
I.color_range = quantile(X, [0.025, 0.975]);
I.colormap = parula(64);
I = parse_optInputs_keyvalue(varargin, I);

G_plot = G;
G_plot.grid_data{1} = flipud(rot90(G.grid_data{1}));
G_plot.grid_data{2} = rot90(G.grid_data{2},3);

% plot p map on the downsampled surface
figure;
subplot(1,2,1);
imagesc(G_plot.grid_data{1}, I.color_range);
colormap(I.colormap);
colorbar;
title('Right Hemi');
subplot(1,2,2);
imagesc(G_plot.grid_data{2}, I.color_range); %#ok<FLUDLR>
colormap(I.colormap);
title('Left Hemi');
colorbar;