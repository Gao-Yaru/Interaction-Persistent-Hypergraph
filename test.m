clear;

% [n,Data] = hypergraph.GetData('baboon');
% RawFiltration = hypergraph.GetFiltration(n,Data);
% Filtration = hypergraph.ModFiltration(RawFiltration);
% [Fil1,SFil1,SFil2] = hypergraph.GetSortedFiltration(Filtration);
% [BM0,BM1] = hypergraph.GetBoundaryMatrix(n,Fil1,SFil1,SFil2);
% P0 = hypergraph.GetPivots(BM0);
% P1 = hypergraph.GetPivots(BM1);
% [bar0,bar1,h1] = hypergraph.GetBarcodes(P0,P1,Fil1,SFil1,SFil2,Filtration);
% 
% K = nchoosek(1:13,3);
% M = zeros(144,3);
% for i = 1:144
%     M(i,:) = K(SFil2(i,1),:);
% end

[bar0,bar1,h1] = hypergraph.Barcodes('baboon');
% [bar0,bar1,h1] = hypergraph.Barcodes('malawi');
% [bar0,bar1,h1] = hypergraph.Barcodes('confer');
% [bar0,bar1,h1] = hypergraph.Barcodes('indust');
% [bar0,bar1,h1] = hypergraph.Barcodes('highsc');
% 
% c = 0;
% [n,m] = size(bar1);
% for i = 1:n
%     if bar1(i,2) > 100
%         c = c + 1;
%     end
% end
% d = 0;
% [n,m] = size(h1);
% for i = 1:n
%     if h1(i,2) > 100
%         d = d + 1;
%     end
% end

% [bar0,bar1,h1] = hypergraph.Barcodes('hospit');

% plot_bars(bar0',0,0,10,0);
% plot_bars(bar1',1,0,10,0);
% plot_bars(h1',1,0,10,1);