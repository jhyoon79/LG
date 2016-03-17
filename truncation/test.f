c234567
        external func
        real*8 x,s
        call qsimp(func(x),1,3,s)
        print *,s
        stop
        end
        
        real*8 function func(x)
        real*8 x
        func = x**2.
        return
        end 
