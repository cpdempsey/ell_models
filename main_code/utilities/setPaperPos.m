function setPaperPos(position_vec)

set(gcf,'position',position_vec,'paperunits','points','paperposition',...
    [0 0 position_vec(3) position_vec(4)],'papersize',[position_vec(3) position_vec(4)])

end