function [res_x, idx_of_result] = knee_pt(y, x, just_return)
    % Returns the x-location of a knee in the curve y=f(x).
    % The knee is determined by fitting two lines, one to the left 
    % and one to the right of a bisection point. The knee is the point
    % which minimizes the combined error of these two fits.
    %
    % Parameters:
    % y          - Input curve (vector with >=3 elements)
    % x          - Corresponding x-values for y (same size as y)
    % just_return- If true, the function will return NaN on errors
    %              without throwing them
    %
    % Returns:
    % res_x       - x-value of the knee point
    % idx_of_result - Index of x where the knee occurs
    %
    % If x is not specified or empty, it defaults to 1:length(y).
    % If just_return is not specified or empty, it defaults to false.
	
	% Dmitry Kaplan (2023). Knee Point (https://www.mathworks.com/matlabcentral/fileexchange/35094-knee-point), 
	% MATLAB Central File Exchange. Retrieved October 30, 2023.
	
	
	
    
    % Internal flag for error calculation method
    use_absolute_dev_p = true; % If false, uses quadratic error

    % Initialize default values and error-handling mode
    res_x = NaN;
    idx_of_result = NaN;
    issue_errors = ~(nargin > 2 && ~isempty(just_return) && just_return);

    % Basic validation of input y
    if isempty(y)
        handleError('y cannot be empty', issue_errors);
        return;
    end
    if sum(size(y) == 1) ~= 1
        handleError('y must be a vector', issue_errors);
        return;
    end
    y = y(:);
    if length(y) < 3
        handleError('y must be at least 3 elements long', issue_errors);
        return;
    end

    % Set x-values if not provided
    if nargin < 2 || isempty(x)
        x = (1:length(y))';
    else
        x = x(:);
    end

    % Ensure x and y have matching dimensions
    if any(size(x) ~= size(y))
        handleError('x and y must have the same dimensions', issue_errors);
        return;
    end

    % Sort data by x-values if necessary
    if any(diff(x) < 0)
        [x, sorting_indices] = sort(x);
        y = y(sorting_indices);
    else
        sorting_indices = 1:length(x);
    end

    % Compute regression parameters for both left and right of the potential knee
    [slope_left, intercept_left] = computeRegressionParams(x, y, 'forward');
    [slope_right, intercept_right] = computeRegressionParams(x, y, 'backward');

    % Compute cumulative error for each potential knee point
    error_curve = computeErrorCurve(x, y, slope_left, intercept_left, slope_right, intercept_right, use_absolute_dev_p);

    % Identify the knee as the point with the minimum cumulative error
    [~, loc] = min(error_curve);
    res_x = x(loc);
    idx_of_result = sorting_indices(loc);
end

function [slope, intercept] = computeRegressionParams(x, y, direction)
    % Compute linear regression parameters for data.
    % 'direction' can be 'forward' (from start to end) or 'backward' (from end to start).
    switch direction
        case 'forward'
            sigma_xy = cumsum(x .* y);
            sigma_x = cumsum(x);
            sigma_y = cumsum(y);
            sigma_xx = cumsum(x.^2);
            n = (1:length(y))';
        case 'backward'
            sigma_xy = cumsum(flip(x) .* flip(y));
            sigma_x = cumsum(flip(x));
            sigma_y = cumsum(flip(y));
            sigma_xx = cumsum(flip(x).^2);
            n = (1:length(y))';
        otherwise
            error('Invalid direction specified');
    end
    determinant = n .* sigma_xx - sigma_x.^2;
    slope = (n .* sigma_xy - sigma_x .* sigma_y) ./ determinant;
    intercept = (sigma_xx .* sigma_y - sigma_x .* sigma_xy) ./ determinant;

    if strcmp(direction, 'backward')
        slope = flip(slope);
        intercept = flip(intercept);
    end
end

function error_curve = computeErrorCurve(x, y, slope_left, intercept_left, slope_right, intercept_right, use_absolute)
    % Compute the error curve based on given slopes and intercepts.
    error_curve = NaN(size(y));
    for breakpt = 2:(length(y) - 1)
        left_errors = (slope_left(breakpt) * x(1:breakpt) + intercept_left(breakpt)) - y(1:breakpt);
        right_errors = (slope_right(breakpt) * x(breakpt:end) + intercept_right(breakpt)) - y(breakpt:end);

        if use_absolute
            error_curve(breakpt) = sum(abs(left_errors)) + sum(abs(right_errors));
        else
            error_curve(breakpt) = sqrt(sum(left_errors.^2)) + sqrt(sum(right_errors.^2));
        end
    end
end

function handleError(msg, issue_errors)
    % Handle errors: either throw them or return gracefully.
    if issue_errors
        error(['knee_pt: ', msg]);
    end
end
