classdef hypergraph
    % compute barcodes of dim 0 and 1 for hypergraphs from data set files
    methods(Static)
        
        % get data from file
        function [n,Data] = GetData(t)
            % high school 2013 file
            if t == 'highsc'
                n = 327;
                fid = fopen('HighSchool2013.csv');
                File = textscan(fid, '%f %f %f %s %s', 'delimiter', ',');
                A = File{1,1};
                B = File{1,2};
                C = File{1,3};
                Pop = zeros(n,1);
                j = 1;
                for i = 1:length(A)
                    if ismember(B(i,1),Pop) == 0
                        Pop(j,1) = B(i,1);
                        j = j + 1;
                    end
                    if ismember(C(i,1),Pop) == 0
                        Pop(j,1) = C(i,1);
                        j = j + 1;
                    end
                end
                Data = zeros(length(A),3);
                for i = 1:length(A)
                    Data(i,1) = A(i,1) - 1385982020;
                    for j = 1:n
                        if B(i,1) == Pop(j,1)
                            Data(i,2) = j;
                        elseif C(i,1) == Pop(j,1)
                            Data(i,3) = j;
                        end
                    end
                end
            end
            
            % workplace 2013 file
            if t == 'indust'
                n = 217;
                fid = fopen('InVS15.txt');
                File = textscan(fid, '%f %f %f', 'delimiter', ' ');
                A = File{1,1};
                B = File{1,2};
                C = File{1,3};
                Data = zeros(length(A),3);
                Pop = zeros(n,1);
                j = 1;
                for i = 1:length(A)
                    if ismember(B(i,1),Pop) == 0
                        Pop(j,1) = B(i,1);
                        j = j + 1;
                    end
                    if ismember(C(i,1),Pop) == 0
                        Pop(j,1) = C(i,1);
                        j = j + 1;
                    end
                end
                for i = 1:length(A)
                    Data(i,1) = A(i,1) - 28820;
                    for j = 1:n
                        if B(i,1) == Pop(j,1)
                            Data(i,2) = j;
                        end
                        if C(i,1) == Pop(j,1)
                            Data(i,3) = j;
                        end
                    end
                end
            end
            
            % hospital file
            if t == 'hospit'
                n = 73;
                fid = fopen('LH10.txt');
                File = textscan(fid, '%f %f %f', 'delimiter', ' ');
                A = File{1,1};
                B = File{1,2};
                C = File{1,3};
                Data = zeros(length(A),3);
                Pop = zeros(n,1);
                j = 1;
                for i = 1:length(A)
                    if ismember(B(i,1),Pop) == 0
                        Pop(j,1) = B(i,1);
                        j = j + 1;
                    end
                    if ismember(C(i,1),Pop) == 0
                        Pop(j,1) = C(i,1);
                        j = j + 1;
                    end
                end
                for i = 1:length(A)
                    Data(i,1) = A(i,1);
                    for j = 1:n
                        if B(i,1) == Pop(j,1)
                            Data(i,2) = j;
                        end
                        if C(i,1) == Pop(j,1)
                            Data(i,3) = j;
                        end
                    end
                end
            end
            
            % conference file
            if t == 'confer'
                n = 403;
                fid = fopen('SFHH.txt');
                File = textscan(fid, '%f %f %f', 'delimiter', ' ');
                A = File{1,1};
                B = File{1,2};
                C = File{1,3};
                Data = zeros(length(A),3);
                Pop = zeros(n,1);
                j = 1;
                for i = 1:length(A)
                    if ismember(B(i,1),Pop) == 0
                        Pop(j,1) = B(i,1);
                        j = j + 1;
                    end
                    if ismember(C(i,1),Pop) == 0
                        Pop(j,1) = C(i,1);
                        j = j + 1;
                    end
                end
                for i = 1:length(A)
                    Data(i,1) = A(i,1) - 32520;
                    for j = 1:n
                        if B(i,1) == Pop(j,1)
                            Data(i,2) = j;
                        end
                        if C(i,1) == Pop(j,1)
                            Data(i,3) = j;
                        end
                    end
                end
            end
            
            % malawi file
            if t == 'malawi'
                n = 86;
                fid = fopen('malawi.csv');
                File = textscan(fid, '%f %f %f %f %f', 'delimiter', ',', 'headerlines',1);
                A = File{1,2};
                B = File{1,4};
                C = File{1,5};
                Data = zeros(length(A),3);
                for i = 1:length(A)
                    Data(i,1) = A(i,1);
                    Data(i,2) = B(i,1);
                    Data(i,3) = C(i,1);
                end
            end
            
            % baboon file
            if t == 'baboon'
                n = 13;
                fid = fopen('RFID_data.txt');
                File = textscan(fid, '%f %s %s %s', 'delimiter', '	', 'headerlines',1);
                A = File{1,1};
                B = File{1,2};
                C = File{1,3};
                
                Pop = {};
                i = 1;
                while (length(Pop) < n) && (i <= length(A))
                    a = 0;
                    for j = 1:length(Pop)
                        if (isequal(B{i,1}, Pop{j,1}) == 1)
                            a = 1;
                            break
                        end
                    end
                    if a == 0
                        Pop(length(Pop) + 1,1) = {B{i,1}};
                    end
                    for j = 1:length(Pop)
                        if (isequal(C{i,1}, Pop{j,1}) == 1)
                            a = 1;
                            break
                        end
                    end
                    if a == 0
                        Pop(length(Pop) + 1,1) = {C{i,1}};
                    end
                    i = i + 1;
                end
                
                Data = zeros(length(A),3);
                for i = 1:length(A)
                    Data(i,1) = A(i,1) - 1560396500;
                    for j = 1:n
                        if (isequal(B{i,1}, Pop{j,1}) == 1)
                            Data(i,2) = j;
                        end
                        if (isequal(C{i,1}, Pop{j,1}) == 1)
                            Data(i,3) = j;
                        end
                    end
                end
            end
        end
        
        % get filtration from data
        function result = GetFiltration(n, A)
            % get filtration from number of vertices n and data A
            
            % initialization and dim 0
            result = {};
            DM = zeros(n,1);
            result(1,1) = {DM};
            
            % dim 1
            DM = zeros(n,n);
            for i = 1:length(A)
                if A(i,2) < A(i,3)
                    DM(A(i,2),A(i,3)) = DM(A(i,2),A(i,3)) + 1;
                else
                    DM(A(i,3),A(i,2)) = DM(A(i,3),A(i,2)) + 1;
                end
            end
            result(2,1) = {DM};
            
            % dim 2
            DM = {};
            for i = 1:n
                DM(i,1) = {zeros(n,n)};
            end
            % at each particular time
            p1 = 1;
            p2 = 1;
            while p2 < length(A) + 1
                p1 = p2;
                a = A(p1,1);
                while (p2 < length(A) + 1) && (A(p2,1) == a)
                    p2 = p2 + 1;
                end
                % record connection matrix
                M = zeros(n,n);
                for i = p1:p2 - 1
                    M(A(i,2),A(i,3)) = M(A(i,2),A(i,3)) + 1;
                    M(A(i,3),A(i,2)) = M(A(i,3),A(i,2)) + 1;
                end
                % record dim 2 complex
                % record removed dim 1 items
                removed = zeros(n,n);
                for i = 1:n
                    C = zeros(n,1);
                    c = 0;
                    for j = i + 1:n
                        if M(i,j) == 1
                            c = c + 1;
                            C(c,1) = j;
                        end
                    end
                    if c > 1
                        for a = 1:c
                            for b = a + 1:c
                                if M(C(a,1),C(b,1)) == 1
                                    DM{i,1}(C(a,1),C(b,1)) = DM{i,1}(C(a,1),C(b,1)) + 1;
                                    if removed(i,C(a,1)) == 0
                                        result{2,1}(i,C(a,1)) = result{2,1}(i,C(a,1)) - 1;
                                        removed(i,C(a,1)) = 1;
                                    end
                                    if removed(i,C(b,1)) == 0
                                        result{2,1}(i,C(b,1)) = result{2,1}(i,C(b,1)) - 1;
                                        removed(i,C(b,1)) = 1;
                                    end
                                    if removed(C(a,1),C(b,1)) == 0
                                        result{2,1}(C(a,1),C(b,1)) = result{2,1}(C(a,1),C(b,1)) - 1;
                                        removed(C(a,1),C(b,1)) = 1;
                                    end
                                end
                            end
                        end
                    end
                end
            end
            result(3,1) = {DM};
        end
        
        % modify the filtrations
        % t -> log(max(t)/t)+1
        function result = ModFiltration(RawFiltration)
            result = {};
            % dim 0
            result(1,1) = {RawFiltration{1,1}};
            % dim 1
            M = RawFiltration{2,1};
            a = max(max(M));
            [n,m] = size(M);
            for i = 1:n
                for j = i + 1:m
                    if M(i,j) == 0
                        M(i,j) = 1e21;
                    else
                        M(i,j) = log(a / M(i,j))+1;
                    end
                end
            end
            result(2,1) = {M};
            % dim 2
            M = RawFiltration{3,1};
%             A = zeros(n,1);
%             for i = 1:n
%                 A(i,1) = max(max(M{i,1}));
%             end
%             a = max(A);
            for i = 1:n
                for j = i + 1:n
                    for k = j + 1:n
                        if M{i,1}(j,k) == 0
                            M{i,1}(j,k) = 1e21;
                        else
                            M{i,1}(j,k) = log(a / M{i,1}(j,k))+1;
                        end
                    end
                end
            end
            result(3,1) = {M};
        end
        
        % sort the filtrations
        % Fil1, sorted decreasingly
        % SFil1: sorted increasingly, remove infinities
        % SFil2: sorted increasingly, remove infinities
        function [Fil1,SFil1,SFil2] = GetSortedFiltration(Filtration)
            n = length(Filtration{1,1});
            % dim 1
            M = zeros(nchoosek(n,2),1);
            c = 0;
            for i = 1:n
                for j = i + 1:n
                    c = c + 1;
                    M(c,1) = Filtration{2,1}(i,j);
                end
            end
            [N,Fil1] = sort(M(:,1),'descend');
            i = 1;
            I = nchoosek(1:n,2);
            while (i <= nchoosek(n,2)) && (Filtration{2,1}(I(Fil1(i,1),1),I(Fil1(i,1),2)) > 1e20)
                i = i + 1;
            end
            SFil1 = flip(Fil1(i:nchoosek(n,2),1));
            % dim 2
            N = zeros(nchoosek(n,3),1);
            c = 0;
            for i = 1:n
                for j = i + 1:n
                    for k = j + 1:n
                        c = c + 1;
                        N(c,1) = Filtration{3,1}{i,1}(j,k);
                    end
                end
            end
            [M,Fil2] = sort(N(:,1),'descend');
            i = 1;
            I = nchoosek(1:n,3);
            while (i <= nchoosek(n,3)) && (Filtration{3,1}{I(Fil2(i,1),1),1}(I(Fil2(i,1),2),I(Fil2(i,1),3)) > 1e20)
                i = i + 1;
            end
            SFil2 = flip(Fil2(i:nchoosek(n,3),1));
        end
        
        % get boundary matrix
        function [BM0,BM1] = GetBoundaryMatrix(n,Fil1,SFil1,SFil2)
            % dim 0
            k1 = length(Fil1);
            K1 = nchoosek(1:n,2);
            BM0 = zeros(n,length(SFil1));
            for i = 1:length(SFil1)
                BM0(K1(SFil1(i,1),1),i) = 1;
                BM0(K1(SFil1(i,1),2),i) = 1;
            end
            % dim 1
            K2 = nchoosek(1:n,3);
            M = zeros(n,n);
            c = 1;
            for i = 1:n
                for j = i + 1:n
                    M(i,j) = c;
                    c = c + 1;
                end
            end
            N = zeros(k1,1);
            for i = 1:k1
                N(Fil1(i,1),1) = i;
            end
            BM1 = zeros(k1,length(SFil2));
            for i = 1:length(SFil2)
                BM1(N(M(K2(SFil2(i,1),1),K2(SFil2(i,1),2)),1),i) = 1;
                BM1(N(M(K2(SFil2(i,1),1),K2(SFil2(i,1),3)),1),i) = 1;
                BM1(N(M(K2(SFil2(i,1),2),K2(SFil2(i,1),3)),1),i) = 1;
            end
        end
        
        % get pivots from boundary matrix
        function P = GetPivots(BM)
            [a,b] = size(BM);
            SBM = BM;
            P = zeros(a,2);
            A = zeros(b,1); % remembers marked columns
            for i = 1:a
                P(i,1) = i;
            end
            for i = 1:a
                j = 1;
                while (j <= b) && ((SBM(i,j) == 0) | (A(j,1) == 1))
                    j = j + 1;
                end
                if j <= b
                    P(i,2) = j;
                    A(j,1) = 1;
                    for k = j + 1:b
                        if SBM(i,k) == 1
                            for m = i:a
                                SBM(m,k) = xor(SBM(m,k),SBM(m,j));
                            end
                        end
                    end
                end
            end
        end
        
        % get barcodes from pivots and filtration
        function [barcode0,barcode1,Hhat1] = GetBarcodes(P0,P1,Fil1,SFil1,SFil2,Filtration)
            % dim 0
            [n,m] = size(P0);
            b0 = length(SFil1);
            barcode0 = zeros(n,2);
            K1 = nchoosek(1:n,2);
            for i = 1:n
                if P0(i,2) == 0
                    barcode0(i,2) = 1e21;
                else
                    barcode0(i,2) = Filtration{2,1}(K1(SFil1(P0(i,2),1),1),K1(SFil1(P0(i,2),1),2));
                end
            end
            
            % dim 1
            barcode1 = zeros(length(Fil1),2);
            Hhat1 = zeros(length(Fil1),2);
            K2 = nchoosek(1:n,3);
            for i = 1:length(Fil1)
                a = Filtration{2,1}(K1(Fil1(i,1),1),K1(Fil1(i,1),2));
                if P1(i,2) == 0
                    b = 1e21;
                else
                    b = Filtration{3,1}{K2(SFil2(P1(i,2),1),1),1}(K2(SFil2(P1(i,2),1),2),K2(SFil2(P1(i,2),1),3));
                end
                if a < b
                    barcode1(i,1) = a;
                    barcode1(i,2) = b;
                else
                    if a > b
                        Hhat1(i,1) = b;
                        Hhat1(i,2) = a;
                    end
                end
            end
            barcode1 = barcode1(any(barcode1,2),:);
            Hhat1 = Hhat1(any(Hhat1,2),:);
        end
        
        % get barcodes from file
        function [bar0,bar1,h1] = Barcodes(t)
            [n,Data] = hypergraph.GetData(t);
            RawFiltration = hypergraph.GetFiltration(n, Data);
            Filtration = hypergraph.ModFiltration(RawFiltration);
            [Fil1,SFil1,SFil2] = hypergraph.GetSortedFiltration(Filtration);
            [BM0,BM1] = hypergraph.GetBoundaryMatrix(n,Fil1,SFil1,SFil2);
            P0 = hypergraph.GetPivots(BM0);
            P1 = hypergraph.GetPivots(BM1);
            [bar0,bar1,h1] = hypergraph.GetBarcodes(P0,P1,Fil1,SFil1,SFil2,Filtration);
        end
    end
end