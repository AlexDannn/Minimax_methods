function [ val ] = r_res( x, a, b )
    top_val = 0;
    bot_val = 0;
        for i = 1:length(a)
            top_val = top_val + a(i)*x^(i-1);
        end
        for j = 1:length(b)
            bot_val = bot_val + b(j)*x^(j-1);
        end
    val = top_val/bot_val;
end

