#include <stdio.h>
#include <stdlib.h>

#include "CommandLines.h"
#include "Process_Read.h"
#include "Assembly.h"
#include "Levenshtein_distance.h"
#include "htab.h"
#include "hifiasm_entry.h"

int hifiasm_main(int argc, char *argv[])
{
    int i, ret;

    yak_reset_realtime();

    init_opt(&asm_opt);

    /*
     * hifiasm 原有参数解析逻辑保持不变。
     * --help、--version 或参数不足时，CommandLine_process()
     * 可能返回 false。
     */
    const int parse_result = CommandLine_process(argc, argv, &asm_opt);
    if (parse_result <= 0)
    {
        /*
         * Keep the original hifiasm main() behaviour here.  Parameter
         * parsing may stop before every owned field has been initialized,
         * so destory_opt() is only safe after an assembly run.
         */
        return parse_result < 0 ? 1 : 0;
    }

    if (asm_opt.sec_in)
    {
        ret = ha_assemble_pair();
    }
    else if (asm_opt.dbg_ovec_cal)
    {
        ret = ha_ec_dbg();
    }
    else
    {
        ret = ha_assemble();
    }

    destory_opt(&asm_opt);

    fprintf(stderr, "[M::%s] Version: %s\n",
            __func__, HA_VERSION);

    fprintf(stderr, "[M::%s] CMD:", __func__);

    for (i = 0; i < argc; ++i)
    {
        fprintf(stderr, " %s", argv[i]);
    }

    fprintf(stderr,
            "\n[M::%s] Real time: %.3f sec; "
            "CPU: %.3f sec; Peak RSS: %.3f GB\n",
            __func__,
            yak_realtime(),
            yak_cputime(),
            yak_peakrss_in_gb());

    return ret;
}

/*
 * 编译原始独立 hifiasm 时定义：
 *
 *     -DHIFIASM_STANDALONE
 *
 * 集成进 HapFold 时不定义该宏，因此最终只有一个 main。
 */
#ifdef HIFIASM_STANDALONE

int main(int argc, char *argv[])
{
    return hifiasm_main(argc, argv);
}

#endif
