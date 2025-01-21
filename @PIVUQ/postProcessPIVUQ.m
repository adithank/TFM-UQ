
function output = postProcessPIVUQ(output)
    for i = 1: length(output)
        
        [X,Y] = meshgrid(output(i).xvec,output(i).yvec);
        output(i).UPost = output(i).U; output(i).VPost = output(i).V;
        output(i).UStd = squeeze(nanstd(output(i).Usamp,0,3)); output(i).VStd = squeeze(nanstd(output(i).Vsamp,0,3));

        delIndsT = [];

        for j = 1 : size(output(i).U,3)
            [xdb, ydb, Udb, Vdb, delInds,USdb,VSdb] = dbscan_clusters_UQ(X,Y,squeeze(output(i).U(:,:,j)),squeeze(output(i).V(:,:,j)),squeeze(output(i).Usamp(:,:,j,:)),squeeze(output(i).Vsamp(:,:,j,:)));
            delIndsT(:,:,j) = (delInds);
        end
        delIndsT = logical(delIndsT);
        output(i).delInds = delIndsT;
        output(i).UPost(delIndsT) = nan; output(i).VPost(delIndsT) = nan;
        output(i).UStd(delIndsT) = nan; output(i).VStd(delIndsT) = nan;
    end
end




function [xdb, ydb, Udb, Vdb, delInds,USdb,VSdb] = dbscan_clusters_UQ(xvec,yvec,U,V,US,VS);
dx = sort(unique(xvec(:)));
dx = dx(2) - dx(1);

% [xvec,yvec] = meshgrid(xvec,yvec);

nB = size(US,3);

xdb = []; ydb = []; Udb = []; Vdb = [];
USdb = []; VSdb = [];
delInds = logical(zeros(size(US,[1 2])));

for i = 1: size(xvec,1)
    for j = 1 : size(yvec,2)

        XX = [squeeze(US(i,j,:)) squeeze(VS(i,j,:))];
        nanCount = sum(sum(isnan(XX),2) ~= 0);

        if nanCount == (nB)
            delInds(i,j) = true;
            xdb(end+1) = xvec(i,j);
            ydb(end+1) = yvec(i,j);
            Udb(end+1) = nan;
            Vdb(end+1) = nan;
            USdb(end+1) = nan;
            VSdb(end+1) = nan;
            continue;
        end

        % Find clusters using dbscan
        idx = dbscan(XX, dx*0.2, ceil(nB/10));

        if (sum(idx==-1)+nanCount) > (0.3 * nB) % Too many outliers
            delInds(i,j) = true;
            xdb(end+1) = xvec(i,j);
            ydb(end+1) = yvec(i,j);
            Udb(end+1) = nan;
            Vdb(end+1) = nan;
            USdb(end+1) = nan;
            VSdb(end+1) = nan;

        elseif length(unique(idx(idx~=-1))) == 1
            xdb(end+1) = xvec(i,j);
            ydb(end+1) = yvec(i,j);
            Udb(end+1) = U(i,j);
            Vdb(end+1) = V(i,j);
            USdb(end+1) = std(US(i,j,:),0,'all');
            VSdb(end+1) = std(VS(i,j,:),0,'all');

        elseif (length(unique(idx(~(idx==-1)))) ~= 1) % If multiple, acceptable clusters
            delInds(i,j) = true;

            ids =  unique(idx(idx~=-1));
            for id = ids'
                xdb(end+1) = xvec(i,j);
                ydb(end+1) = yvec(i,j);
                Udb(end+1) = mean(US(i,j,idx==id),'all');
                Vdb(end+1) = mean(VS(i,j,idx==id),'all');
                USdb(end+1) = std(US(i,j,idx==id),0,'all');
                VSdb(end+1) = std(VS(i,j,idx==id),0,'all');
            end
        end
    end
end
end
