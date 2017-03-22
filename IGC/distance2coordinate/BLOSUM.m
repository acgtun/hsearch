function BLOSUM
    % converted distance matrix
    D = [0 11 14 14 13 11 11 10 16 10 10 11 11 14 13 6 9 21 15 8 
         11 0 11 15 20 8 10 15 13 15 13 6 12 17 16 11 12 22 16 15 
         14 11 0 10 21 11 11 12 12 16 16 11 15 18 17 8 11 25 17 16 
         14 15 10 0 21 11 7 14 16 16 18 13 17 18 15 10 13 25 19 16 
        13 20 21 21 0 20 22 21 23 15 15 20 16 19 22 15 16 24 20 15 
        11 8 11 11 20 0 6 15 13 15 13 8 10 17 14 9 12 20 14 13 
        11 10 11 7 22 6 0 15 13 15 15 8 14 17 14 9 12 22 16 13 
        10 15 12 14 21 15 15 0 18 18 18 15 17 18 17 10 15 21 19 16 
        16 13 12 16 23 13 13 18 0 18 18 15 17 16 19 14 17 23 11 18 
        10 15 16 16 15 15 15 18 18 0 4 15 7 10 17 12 11 21 13 2 
        10 13 16 18 15 13 15 18 18 4 0 13 5 10 17 12 11 19 13 6 
        11 6 11 13 20 8 8 15 15 15 13 0 12 17 14 9 12 22 16 13 
        11 12 15 17 16 10 14 17 17 7 5 12 0 11 16 11 12 18 14 7 
        14 17 18 18 19 17 17 18 16 10 10 17 11 0 21 14 15 15 7 12 
        13 16 17 15 22 14 14 17 19 17 17 14 16 21 0 13 14 26 20 15 
        6 11 8 10 15 9 9 10 14 12 12 9 11 14 13 0 7 21 15 12 
        9 12 11 13 16 12 12 15 17 11 11 12 12 15 14 7 0 20 16 9 
        21 22 25 25 24 20 22 21 23 21 19 22 18 15 26 21 20 0 14 21 
        15 16 17 19 20 14 16 19 11 13 13 16 14 7 20 15 16 14 0 13 
        8 15 16 16 15 13 13 16 18 2 6 13 7 12 15 12 9 21 13 0]

	for d = 1 : 50
        Y = mdscale(D, d);
        pd = DistanceMatrix(Y);
        e1(d) = L1error(D, pd);
        e2(d) = L2error(D, pd);
        fprintf('d = %d e1 = %f e2 = %f\n', d, e1(d), e2(d));
    end
    plot([1 : 50], e1, 'r+-')
    hold on
    plot([1 : 50], e2, 'b*-')
end

function e = L2error(D, pd)
    e = 0;
    for i = 1 : 20
        for j = 1 : 20
            e = e + (D(i, j) - pd(i, j))* (D(i, j) - pd(i, j));
        end
    end
end

function e = L1error(D, pd)
    e = 0;
    for i = 1 : 20
        for j = 1 : 20
            e = e + abs(D(i, j) - pd(i, j));
        end
    end
end

function pd = DistanceMatrix(Y)
    for i = 1 : 20
        for j = 1 : 20
            a = Y(i, : );
            b = Y(j, : );
            pd(i, j) = pdist2(a, b);
        end
    end
end

function PrintConvertedCoordinates(Y)
    MarkColor = [237 125 49]./ 255;
    AA = {'Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val'}
    for i = 1 : 20
        plot3(Y(i, 1), Y(i, 2), Y(i, 3), '*', 'MarkerSize',10, 'MarkerFaceColor', 'r',  'MarkerEdgeColor', 'r');
        hold on
        text(Y(i, 1) + 0.5, Y(i, 2) + 0.5, Y(i, 3) + 0.5 ,AA{i}, 'FontSize', 10);
        hold on
    end
    
    grid on
    set(gca, 'FontSize',20)
end
