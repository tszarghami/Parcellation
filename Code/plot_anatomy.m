%% Plot MIPs of partitions on glass brains

% Code partially adapted from spm_mb_anatomy.m in SPM12

VOX = MB.VOX;
XYZ = MB.XYZ;

nz  = numel(MB.z);                             % number of particles
ns  = size(MB.W,2);                            % number of states (n - 1)
nx  = size(MB.J,2);                            % number of states (n)
col = spm_MB_col(nz);                           % colours

% functional anatomy of states
%--------------------------------------------------------------------------
figure
c    = var(abs(MB.Y));
spm_mip_mod(logical(abs(MB.V))*c',XYZ,VOX);
axis image, axis off
str{1} = sprintf('%d states, %d particles',ns,nz);
str{2} = sprintf('(%d eigenstates)',nx);
title(str,'FontSize',12)
colormap gray

%% Visualize some/all of the partitions from requested scales

for s = scales_to_plot
    IDX_scale = IDX_cell{s};
    N_partitions = size(IDX_scale,1);
    % N_partitions_divis = sum(max(IDX_scale,[],2)==2);
    partition_num = 1;
    for partition = 1:N_partitions
        idx = IDX_scale(partition,:);
        clear ax idx1 c0 c
         % plot both divisible and indivisible partitions
         % but for scale analysis, indivisible partitions belong to previous scale
        if (partition_num < max_partitions)
            h=figure;
            for j=1:max(idx) % sub-partition 1 or 2
                idx1 = (idx==j);
                Y=MB.Y;Y=Y(:,idx1);
                c    = var(abs(Y));
                V=MB.V;V=V(:,idx1);
                ax = subplot(1,max(idx),j);
                spm_mip_mod(logical(abs(V))*c',XYZ,VOX);
                axis image, axis off,
                if any(dot(idx1,IDX_cell{2}==1)) % adjust the colormap based on the first bi-partition
                    colormap (ax,'bone')
                elseif any(dot(idx1,IDX_cell{2}==2))
                    colormap(ax,'pink')
                end
            
                if s>1, title(['Subpartition ',num2str(j)]), end
            end
            t = ['Scale',num2str(s),'-Partition',num2str(partition_num)];
            sgtitle(t) % for github code
            % saveas(h,[t,'.fig'])
            saveas(h,[t,'.png'])
            partition_num = partition_num + 1;
        end
    end
end

%% Hierarchical depth of regions

VOX = MB.VOX;
XYZ = MB.XYZ;
V = MB.V;
IDX_mat = cell2mat(IDX_cell(:));
Depth = sum(IDX_mat>0)-2;
mip = spm_mip_mod(logical(abs(V))*Depth',XYZ,VOX);
mip = (mip/max(mip,[],'all'));
mip = 1-mip; % To cancel out 1-mip in spm_mip (line 114)
mip = mip* max(Depth);

figure
image(mip,'CDataMapping','scaled')
axis image, axis off, 
c = colormap('jet');
c(1,:) = [0.9 0.9 0.9];
colormap(c)
title('Depth of Regions in the Hierarchy')

clb = colorbar(gca,'Ticks',[min(Depth):1:max(Depth)]);
title(clb,'Depth','fontweight','bold')
caxis([min(Depth),max(Depth)])
set(gca,'fontweight','bold','fontsize',14)
saveas(gcf,'Depth.png')

%% Visualize the final parcellation 
% (leaf nodes + parcels with N=1 from prev scales)

IDX_partitions = zeros(1,size(IDX_mat,2));
V = MB.V;
s =1;
for k = 1:numel(N_opt)
    n_opt = N_opt{k};
    partition = (n_opt == 1);
    partition = reshape(partition,[],1);
    idx = IDX_raw{k};
    idx = reshape(idx,[],size(idx,3));
    idx = idx(partition,:);
    for j = 1:size(idx,1)
        IDX_partitions = IDX_partitions+ s*idx(j,:);
        s = s+20;
    end
end

figure
mip = spm_mip_mod(logical(abs(V))*IDX_partitions',XYZ,VOX);
image(mip);
axis image, axis off, 
c = colormap('colorcube');
c(end-30:end,:) = repmat([1,1,1],31,1);
colormap(c);
title('Final Parcellation')
ax = gca;ax.FontWeight = 'bold';ax.FontSize=14;
saveas(gcf,'Parcellation.png')

