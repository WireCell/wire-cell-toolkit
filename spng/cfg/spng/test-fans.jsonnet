local fans = import "spng/fans.jsonnet";

{
    forfilter: fans.fanout_select("test", 3, [["w","g","d"], ["w","g","d"], ["w","g"]])
}
