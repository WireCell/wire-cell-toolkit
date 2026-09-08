// Externally determined "bad" channels.
[
    # CE group: Inactive FE
    4411,  # femb515x12
    4412,  # femb515x13
    9990,  # femb605x10
    11842,  # femb120x03
    # CE group: Broken connection
    1,  # femb311u39
    400,  # femb301u40
    401,  # femb301u39
    800,  # femb320v01
    801,  # femb320v02
    876,  # femb319v37
    1200,  # femb310v01
    2961,  # femb501u39
    5321,  # femb216u39
    5363,  # femb217u37
    6132,  # femb215v13
    7058,  # femb213x03
    7295,  # femb202x01
    7681,  # femb611u39
    8080,  # femb601u40
    8328,  # femb607u32
    8480,  # femb620v01
    9282,  # femb620x03
    9283,  # femb620x04
    9736,  # femb611x25
    9854,  # femb602x02
    10800,  # femb105u40
    11024,  # femb110u16
    11457,  # femb110v18
    11459,  # femb110v20
    11463,  # femb110v24
    11469,  # femb110v30: bad in 4875-185-1500, ok in 5803-76-3200
    11517,  # femb109v38
    11669,  # femb105v30
    11679,  # femb105v40
    12756,  # femb110x44
    12801,  # femb411u39
    13001,  # femb416u39
    13081,  # femb418u39
    # CE group: ENC > 2000e
    4410,  # femb515x11: High noise, no signal 5008-76
    #-----
    # CE group excessive sticky
    #femb318x
    1719,  # femb318x24
    5125,  # femb211u35
    7551,  # femb208x33
    7190,  # femb211x39
    7194,  # femb211x43
    7918,  # femb616u02, sticky pedestal (three peaks)
    #-----
    # CE group: good.
    # femb311
    2,  # femb311u38, no signal
    4,  # femb311u36, very sticky pedestal 5308-76
    1632,  # femb320x33, very sticky pedestal 5308-76
    2169,  # femb302x07, Mostly stuck on one bad code, 5308-76
    2450,  # femb308x14, Very noisy (1000 ADC) in run 5759 (20nov2019)
    3541,  # femb516v22, very sticky--signal near zero half the time (5308-81)
    3543,  # femb516v24, very sticky--signal near zero half the time (5308-81)
    3661,  # femb513v22, most signal near zero (5308-81)
    3663,  # femb513v24, most signal near zero (5308-81)
    4061,  # femb503v22, most signal near zero (5308-81)
    4063,  # femb503v24, most signal near zero (5308-81)
    4141,  # femb501v22, signal near zero half the time (5308-81)
    4143,  # femb501v24, signal sometimes near zero (5308-81)
    4377,  # femb516x26, very sticky pedestal
    4379,  # femb516x28, very sticky pedestal
    4381,  # femb516x30, very sticky pedestal
    4383,  # femb516x32, very sticky pedestal
    4385,  # femb516x34, very sticky pedestal
    4387,  # femb516x36, very sticky pedestal
    4521,  # femb513x26, very sticky pedestal
    4523,  # femb513x28, very sticky pedestal
    4525,  # femb513x30, very sticky pedestal
    4527,  # femb513x32, very sticky pedestal
    4529,  # femb513x34, very sticky pedestal
    4531,  # femb513x36, very sticky pedestal
    4652,  # femb501x36, very sticky pedestal
    4654,  # femb501x34, very sticky pedestal
    4656,  # femb501x32, very sticky pedestal
    4658,  # femb501x30, very sticky pedestal
    4660,  # femb501x28, very sticky pedestal
    4658,  # femb501x26, very sticky pedestal
    4748,  # femb503x36, very sticky pedestal
    4750,  # femb503x34, very sticky pedestal
    4752,  # femb503x32, very sticky pedestal
    4754,  # femb503x30, very sticky pedestal
    4756,  # femb503x28, very sticky pedestal
    4758,  # femb503x26, very sticky pedestal
    5361,  # femb217u39, no signal
    7680,  # femb611u40: No signal in 5308-76, end wire
    8501,  # femb620v22, very sticky pedestal
    8503,  # femb620v24, very sticky pedestal
    8821,  # femb612v22, very sticky pedestal
    8823,  # femb612v24, very sticky pedestal
    9261,  # femb601v22, very sticky pedestal
    9263,  # femb601v24, very sticky pedestal
    9305,  # femb620x26, very sticky pedestal
    9307,  # femb620x28, very sticky pedestal
    9309,  # femb620x30, very sticky pedestal
    9311,  # femb620x32, very sticky pedestal
    9313,  # femb620x34, very sticky pedestal
    9315,  # femb620x36, very sticky pedestal
    9689,  # femb612x26, very sticky pedestal
    9691,  # femb612x28, very sticky pedestal
    9693,  # femb612x30, very sticky pedestal
    9695,  # femb612x32, very sticky pedestal
    9697,  # femb612x34, very sticky pedestal
    9699,  # femb612x36, very sticky pedestal
    9772,  # femb601x26, very sticky pedestal
    9774,  # femb601x28, very sticky pedestal
    9776,  # femb601x30, very sticky pedestal
    9778,  # femb601x32, very sticky pedestal
    9780,  # femb601x34, very sticky pedestal
    9782,  # femb601x36, very sticky pedestal
    10102,  # femb608x42, mostly stuck on one code
    10189,  # femb609x03, mostly stuck on one code
    10697,  # femb102u23, mostly stuck on a few classic codes
    10907,  # femb107u13, mostly stuck on one code
    11203,  # femb116v04, stuck on many classic codes
    11270,  # femb115v31, stuck on many classic codes
    11902,  # femb119x15, stuck on two classic codes
    12324,  # femb101x44, stuck on many classic codes
    12333,  # femb101x35, stuck on many classic codes
    12744,  # femb109x08, stuck on many classic codes
    13363,  # femb405u37, very noisy, nosignal 5308-76-4800
    #-----
    # These 16 channels are an intermitently bad ASIC.
    # Matt W. 19oct2018.
    # femb316u
    200,   # femb316u40
    202,   # femb316u38
    204,   # femb316u36
    206,   # femb316u34
    208,   # femb316u32
    # femb316v
    991,   # femb316v32
    993,   # femb316v34
    995,   # femb316v36
    997,   # femb316v38
    999,   # femb316v40
    # femb316x
    1829,   # femb316x38
    1831,   # femb316x40
    1833,   # femb316x42
    1835,   # femb316x44
    1837,   # femb316x46
    1839    # femb316x48
        #-----
]
