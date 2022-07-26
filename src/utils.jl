# julia 1.6
#
# diverse usful functions 

function ncbfind(A, l, a)
	"""
	find position of a in ordered list A
	borrowed from the TSSOS package
	https://github.com/wangjie212/TSSOS
	"""
    low = 1
    high = l
    while low <= high
        mid = Int(ceil(1/2*(low+high)))
        if A[mid] == a
           return mid
        elseif A[mid] < a
            high = mid-1
        else
            low = mid+1
        end
    end
    return 0
end

ncbfind(A, a) = ncbfind(A, length(A), a)