function h = imagescc(img)
    if iscomplex(img)
        img = abs(img);
        warning('cast complex to abs');
    end
    h = imagesc(rot90(squeeze(img),1));
    axis image;
    colorbar;
    xticks([]); yticks([]);
end