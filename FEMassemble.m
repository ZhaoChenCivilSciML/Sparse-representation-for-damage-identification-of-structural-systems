function DOFs = FEMassemble(ele)
if ele < 7
    switch mod(ele,5)
        case 0 
            node1 = 1;
            node2 = 3;
        case 1
            node1 = 1+2*floor(ele/5);
            node2 = 2+2*floor(ele/5);
        case 2
            node1 = 2;
            node2 = 4;
        case 3
            node1 = 2;
            node2 = 3;
        case 4
            node1 = 1;
            node2 = 4;
    end
elseif ele < 11 && ele > 6
    switch mod(ele,5)        
        case 0
            node1 = 3;
            node2 = 6;
        case 2
            node1 = 4;
            node2 = 5;
        case 3
            node1 = 4;
            node2 = 6;
        case 4
            node1 = 3;
            node2 = 5;
    end
else
    switch mod(ele,5)
        case 0
            node1 = 6+2*(ele/5-3);
            node2 = 8+2*(ele/5-3);
        case 1
            node1 = 6+2*(floor(ele/5)-2);
            node2 = 5+2*(floor(ele/5)-2);
        case 2
            node1 = 5+2*(floor(ele/5)-2);
            node2 = 7+2*(floor(ele/5)-2);
        case 3
            node1 = 5+2*(floor(ele/5)-2);
            node2 = 8+2*(floor(ele/5)-2);
        case 4
            node1 = 6+2*(floor(ele/5)-2);
            node2 = 7+2*(floor(ele/5)-2);
    end
end
DOFs = [2*(node1-1)+1,2*node1,2*(node2-1)+1,2*node2];
