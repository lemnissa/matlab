function [Q,D] = ElectroStaticDipoles(XYZ,R,F)

N = size(XYZ, 2);
B11 = zeros(4*N, 4*N); %matrix of potential coefficients
F = [F; zeros(3 * N, 1)]; % matrix of required potentials and electric fields 

% checks that balls do not overlap
for ii = 1 : N
    for jj = ii + 1 : N
        if norm(XYZ(:,ii) - XYZ(:, jj)) < R(ii) + R(jj)
            error('Balls overlap!')
        end
    end
end

%builds blocks for matrix 4N x 4N

%%% first row, first block (block is matrix NxN)
% same as main matrix in ElectroStaticBalls
for ii = 1 : N
    for jj = ii : N
        if jj == ii
            B_11(ii, ii) = 1/R(ii);
        else
            B_11(ii, jj) = 1/norm(XYZ(:,ii) - XYZ(:,jj));
            B_11(jj, ii) = B_11(ii, jj);
        end
    end
end

%%% first row, second block
% \varphi = dot(d, r)/abs(r)^3
for ii = 1 : N
    for jj = ii : N
        if jj == ii
            B_12(ii, ii) = 0;
        else
            B_12(ii, jj) = (XYZ(1, jj) - XYZ(1, ii))/norm(XYZ(:,jj) - XYZ(:,ii))^3;
            B_12(jj, ii) = -B_12(ii, jj); %equation depends on vector of distance 
        end                               %between balls, but r_12 = -r_21, so
    end                                   %somewhere you'll see minuses
end

%%% first row, third block 
for ii = 1 : N
    for jj = ii : N
        if jj == ii
            B_13(ii, ii) = 0;
        else
            B_13(ii, jj) = (XYZ(2, jj) - XYZ(2, ii))/norm(XYZ(:,jj) - XYZ(:,ii))^3;
            B_13(jj, ii) = -B_13(ii, jj);
        end
    end
end

%%% first row, fourth block 
for ii = 1 : N
    for jj = ii : N
        if jj == ii
            B_14(ii, ii) = 0;
        else
            B_14(ii, jj) = (XYZ(3, jj) - XYZ(3, ii))/norm(XYZ(:,jj) - XYZ(:,ii))^3;
            B_14(jj, ii) = -B_14(ii, jj);
        end
    end
end

%%% second row, first block
% E = q*r/abs(r)^3
for ii = 1 : N
    for jj = ii : N
        if jj == ii
            B_21(ii, ii) = 0;
        else
            B_21(ii, jj) = (XYZ(1, jj) - XYZ(1, ii))/norm(XYZ(:,jj) - XYZ(:,ii))^3;
            B_21(jj, ii) = -B_21(ii, jj);
        end
    end
end

%%% second row, second block
% field of a dipole of the ball itself is E = d/R^3,
% field of a dipole of the other ball is E = (d, r)r/abs(r)^5 - d/abs(r)^3
for ii = 1 : N
    for jj = ii : N
        if jj == ii
            B_22(ii, ii) = 1/R(ii)^3;
        else
            B_22(ii, jj) = 3/norm(XYZ(:,jj) - XYZ(:,ii))^3 - 1/norm(XYZ(:,jj) - XYZ(:,ii))^5;
            B_22(jj, ii) = B_22(ii, jj);
        end
    end
end

%%% second row, third block
for ii = 1 : N
    for jj = ii : N
        if jj == ii
            B_23(ii, ii) = 0;
        else
            B_23(ii, jj) = 3*(XYZ(2, jj) - XYZ(2, ii))/norm(XYZ(:,jj) - XYZ(:,ii))^3;
            B_23(jj, ii) = -B_23(ii, jj);
        end
    end
end

%%% second row, fourth block
for ii = 1 : N
    for jj = ii : N
        if jj == ii
            B_24(ii, ii) = 0;
        else
            B_24(ii, jj) = 3*(XYZ(3, jj) - XYZ(3, ii))/norm(XYZ(:,jj) - XYZ(:,ii))^3;
            B_24(jj, ii) = -B_24(ii, jj);
        end
    end
end

%%% third row, first block
for ii = 1 : N
    for jj = ii : N
        if jj == ii
            B_31(ii, ii) = 0;
        else
            B_31(ii, jj) = (XYZ(2, jj) - XYZ(2, ii))/norm(XYZ(:,jj) - XYZ(:,ii))^3;
            B_31(jj, ii) = -B_31(ii, jj);
        end
    end
end

%%% third row, second block
for ii = 1 : N
    for jj = ii : N
        if jj == ii
            B_32(ii, ii) = 0;
        else
            B_32(ii, jj) = 3*(XYZ(1, jj) - XYZ(1, ii))/norm(XYZ(:,jj) - XYZ(:,ii))^3;
            B_32(jj, ii) = -B_32(ii, jj);
        end
    end
end

%%% third row, third block
for ii = 1 : N
    for jj = ii : N
        if jj == ii
            B_33(ii, ii) = 1/R(ii)^3;
        else
            B_33(ii, jj) = 3/norm(XYZ(:,jj) - XYZ(:,ii))^3 - 1/norm(XYZ(:,jj) - XYZ(:,ii))^5;
            B_33(jj, ii) = B_33(ii, jj);
        end
    end
end

%%% third row, fourth block
for ii = 1 : N
    for jj = ii : N
        if jj == ii
            B_34(ii, ii) = 0;
        else
            B_34(ii, jj) = 3*(XYZ(3, jj) - XYZ(3, ii))/norm(XYZ(:,jj) - XYZ(:,ii))^3;
            B_34(jj, ii) = -B_34(ii, jj);
        end
    end
end

%%% fourth row, first block
for ii = 1 : N
    for jj = ii : N
        if jj == ii
            B_41(ii, ii) = 0;
        else
            B_41(ii, jj) = (XYZ(3, jj) - XYZ(3, ii))/norm(XYZ(:,jj) - XYZ(:,ii))^3;
            B_41(jj, ii) = -B_41(ii, jj);
        end
    end
end

%%% fourth row, second block
for ii = 1 : N
    for jj = ii : N
        if jj == ii
            B_42(ii, ii) = 0;
        else
            B_42(ii, jj) = (XYZ(1, jj) - XYZ(1, ii))/norm(XYZ(:,jj) - XYZ(:,ii))^3;
            B_42(jj, ii) = -B_42(ii, jj);
        end
    end
end

%%% fourth row, third block
for ii = 1 : N
    for jj = ii : N
        if jj == ii
            B_43(ii, ii) = 0;
        else
            B_43(ii, jj) = (XYZ(2, jj) - XYZ(2, ii))/norm(XYZ(:,jj) - XYZ(:,ii))^3;
            B_43(jj, ii) = -B_43(ii, jj);
        end
    end
end

%%% fourth row, fourth block
for ii = 1 : N
    for jj = ii : N
        if jj == ii
            B_44(ii, ii) = 1/R(ii)^3;
        else
            B_44(ii, jj) = 3/norm(XYZ(:,jj) - XYZ(:,ii))^3 - 1/norm(XYZ(:,jj) - XYZ(:,ii))^5;
            B_44(jj, ii) = B_44(ii, jj);
        end
    end
end

% joins all 16 blocks together 

C = [B_11 B_12 B_13 B_14; B_21 B_22 B_33 B_44; B_31 B_32 B_33 B_44; B_41 B_42 B_43 B_44];


QD = C \ F;
Q = QD(1 : N, 1);
Dx = QD(N + 1 : 2*N, 1)';
Dy = QD(2 * N + 1 : 3*N, 1)';
Dz = QD(3 * N + 1 : 4*N, 1)';
D = -[Dx; Dy; Dz]';

end
