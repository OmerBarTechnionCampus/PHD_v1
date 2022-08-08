function [idx] = SelectingMaxima(VotingMatrix,selectionType)
%SelectingMaxima is a method to choose possible indexes for maxmimum
%probability of concenzus
% 
% selectionType options - 'pareto', 'max', 'max2', '2D_3Sigma','2D_2Sigma', '2DRMS'
%
% output - indexes at the vooting matrix
%%
selectionType = lower(selectionType);
vh3 = sort(unique(VotingMatrix));

switch (selectionType)
    case 'pareto'
        % "counting" - as a pareto signature
        vh3p = [vh3,ones(size(vh3))-2,zeros(size(vh3))];
        for t = 1:size(vh3p,1)
            vh3p(t,2) = numel(find(VotingMatrix == vh3p(t,1)))/numel(VotingMatrix);
            vh3p(t,3) = sum(vh3p(t:-1:1,2));
        end
        %       figure(2), plot(vh3p(:,1),vh3p(:,3)); %pareto
        %       figure(2), hold on;
        %       figure(2), plot(vh3p(:,1),vh3p(:,2)); %precentage
        %       figure(2), title('pareto')
        
        % getting possible options as consenzus of 75% of the vector sample
        %     [am,vm,lm] = ind2sub(size(H3), find(H3>=min( vh3(vh3>0.75*numel(vs)) )) );
        
        th_pareto = min(vh3(vh3p(:,3)>0.971));   % 3d 3-sigma confidence level (Adjustment computations, Error elipses chapeter, Table 19.4, measure of three dimentional position uncertaties)
        idx = find(  VotingMatrix >= th_pareto );
    case 'max'
%         fi = find(  H3>=min( vh3(vh3>vh3(floor(end-numel(vh3)/4)) ) ) );
        idx = find(  VotingMatrix == max(vh3(1:end)) );
    case 'max2' % first 2 max values
        fi = find(  VotingMatrix == max(vh3(1:end)) );
        tempv = vh3(vh3~=max(vh3));
        idx = [ fi; find(  VotingMatrix == max(tempv(1:end)) )];
    %        
    % 2d sigma confidence level (Adjustment computations, Error elipses chapeter, Table 19.3, measure of two dimentional position uncertaties)
    case '2d_3sigma'     
        idx = find(  VotingMatrix >= 0.989 );
    case '2drms'
        idx = find(  VotingMatrix >= 0.982 );
    case '2d_2sigma'
        idx = find(  VotingMatrix >= 0.865 );
    case '2d_95p'
        idx = find(  VotingMatrix >= 0.95 );
end %switch

if (isempty(idx) && ~strcmpi(selectionType,'2d_2sigma'))
    warning('warning in SelectingMaxima :: selection failed trying "2d_2sigma" option');
    idx = SelectingMaxima(VotingMatrix,'2d_2sigma');
    if (isempty(idx) && ~strcmpi(selectionType,'max2'))
        idx = SelectingMaxima(VotingMatrix,'max2');
        warning('warning in SelectingMaxima :: selection failed trying "max2" option');
    end
end

end %SelectingMaxima


%{
% % % % % %     [cmap] = fixColorMap2Items(colormap('jet'),size(vh3p,1));
% % % % % %     [cmap] = reorderColorMap2Items(H3,vh3p(:,1),cmap);
% % % % % %     figure(p*10), scatter3(AZs(ai).*180./pi,VFs(vi),LDs(li),H3(:),cmap);
    % % getting the maxmimum consenzus
    % [am,vm,lm]=ind2sub(size(H3),find(H3==vh3(end)));

    % getting possible options as consenzus of 75% of the vector sample
%     [am,vm,lm] = ind2sub(size(H3), find(H3>=min( vh3(vh3>0.75*numel(vs)) )) );
        fi = find(  H3>=min( vh3(vh3>vh3(floor(end-numel(vh3)/4)) ) ) );
        fi = find(  H3>=min( vh3(vh3>max(vh3)*0.6)));
        th_pareto = min(vh3(vh3p(:,3)>99.5));
        fi = find(  H3>= th_pareto );
        
%}
