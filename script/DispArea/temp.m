t = reshape(raw, 6, 5);
t = repmat(t, 2, 1);
tt = repmat(t(7, :), 6, 1);
t(7:12, :) = tt ./ t(7:12, :);
