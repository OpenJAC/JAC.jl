function main(d)
    bad = 0;  n = 0
    for f in sort(readdir(d))
        endswith(f, ".canon") || continue
        n += 1
        ex = Meta.parseall(read(joinpath(d, f), String); filename=f)
        for a in ex.args
            if a isa Expr && a.head === :error
                println("   PARSE ERROR  ", f, ": ", a);  bad += 1;  break
            end
        end
    end
    println("parsed $n canonical files, $bad with a parse error")
end
main(ARGS[1])
