function trapz(values, dx)

    (sum(values) - values[1]/2 - values[end]/2)*dx

end

function trapz(values, dx, dy)

    temp = zeros(size(values, 1))
    for i = 1:size(values, 1)
      temp[i] = trapz(values[i, :], dy)
    end
    trapz(temp, dx)

end

function trapz(values, dx, dy, dz)

    temp = zeros(size(values, 1));
    for i = 1:size(values, 1)
        temp[i] = trapz(values[i, :, :], dy, dz)
    end
    trapz(temp, dx)

end

function trapz(values, dx, dy, dz, dw)

    temp = zeros(size(values, 1))
    for i = 1:size(values, 1)
        temp[i] = trapz(values[i, :, :, :], dy, dz, dw)
    end
    trapz(temp, dx)

end

function trapz(values, dx, dy, dz, dw, dv)

    temp = zeros(size(values, 1))
    for i = 1:size(values, 1)
      temp[i] = trapz(values[i, :, :, :, :], dy, dz, dw, dv)
    end
    trapz(temp, dx)

end
