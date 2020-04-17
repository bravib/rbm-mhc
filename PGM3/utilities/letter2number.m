%{
function x=letter2number(a)
    switch(a)
        case 'A'
            x=1;
        case 'R'
            x=2;
        case 'N'
            x=3;
        case 'D'
            x=4;
        case 'C' 
            x=5;
        case 'Q'
            x=6;
        case 'E' 
            x=7;
        case 'G'
            x=8;
        case 'H' 
            x=9;
        case 'I'
            x=10;
        case 'L' 
            x=11;
        case 'K'
            x=12;
        case 'M' 
            x=13;
        case 'F'
            x=14;
        case 'P'
            x=15;
        case 'S'
            x=16;
        case 'T' 
            x=17;
        case 'W'
            x=18;
        case 'Y'
            x=19;
        case 'V'
            x=20;
        case '-'
            x=21;
        otherwise
            x=21;
    end
end
%}

% Implementing Jerome convention:
function x=letter2number(a)
    switch(a)
        case 'A'
            x=1;
        case 'C'
            x=2;
        case 'D'
            x=3;
        case 'E'
            x=4;
        case 'F' 
            x=5;
        case 'G'
            x=6;
        case 'H' 
            x=7;
        case 'I'
            x=8;
        case 'K' 
            x=9;
        case 'L'
            x=10;
        case 'M' 
            x=11;
        case 'N'
            x=12;
        case 'P' 
            x=13;
        case 'Q'
            x=14;
        case 'R'
            x=15;
        case 'S'
            x=16;
        case 'T' 
            x=17;
        case 'V'
            x=18;
        case 'W'
            x=19;
        case 'Y'
            x=20;
        case '-'
            x=21;
        otherwise
            x=21;
    end
end

