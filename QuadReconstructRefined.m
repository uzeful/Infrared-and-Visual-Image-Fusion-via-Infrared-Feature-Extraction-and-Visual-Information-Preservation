function bgImg = QuadReconstructRefined(S, img, minDim)
% This function is used to reconstruct the background image using bezier
% interpolation interpolation on the quadtree structure S.
% Implemented by Zhang Yu (uzeful@163.com).

    % double image 
    img = double(img);

    % Matrix M
    M = [
        1   0   0  0;
       -3   3   0  0;
        3  -6   3  0;
       -1   3  -3  1
        ];

    MT = M'; % transform of matrix M

    % preprocess S
    newS = S + (S > 0);

    % max level
    newS = full(newS);
    maxDim = max(newS(:));
    dim = maxDim;

    % pad image
    newS = padarray(newS, [1 1], 'replicate', 'post');
    newImg = padarray(img, [1 1], 'replicate', 'post');

    % temporal reconstructed image
    tempReconstImg = zeros(size(newImg));

    % Begin loop
    while (dim >= minDim + 1)

        len = length(find(newS == dim));
        if len ~= 0
            % Extrat the corresponding blocks at the current level
            [blks, Sind] = qtgetblk(newImg, newS, dim);

            % pixel locations
            subDim = (dim - 1) / 4;
            row = [1, subDim * 2, subDim * 3, subDim * 4 + 1];

            xx = [row; row; row; row];
            yy = xx';
            inds = sub2ind([dim, dim], yy, xx);

            % generate u and v
            u=linspace(0,1,dim); u = u(:); v = u;
            U2 = [ones(dim, 1), u, u .^ 2, u .^3];
            VT2 = [ones(dim, 1), v, v .^ 2, v .^3]';

            % reconstructed blocks
            reBlkSeq = zeros(dim, dim, len);

            for ii = 1 : len
                blockVal = blks( : , : , ii);
                blockVal = blockVal(inds);
                reblkVal = U2 * M * blockVal * MT * VT2;
                reBlkSeq(:,:,ii) = reblkVal;
            end

            % Set the fusion tag as the image indices of the focused blocks
            tempReconstImg = qtsetblk(tempReconstImg,newS,dim,reBlkSeq);
        end
        dim = (dim - 1) / 2 + 1;
    end
    % Final reconstrcucted background image
    bgImg = tempReconstImg(1 : end - 1, 1 : end - 1);
end