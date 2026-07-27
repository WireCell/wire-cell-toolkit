// test_wrap.jsonnet
// Test cases for wrap.jsonnet functionality

local pg = import "pgraph.jsonnet";
local wrapper = import "wrap.jsonnet";

// Helper to create simple test inodes
local make_inode(type, name) = {
    type: type,
    name: name,
    data: {}
};

// Test 1: No wrapping when finder doesn't match type
local test1 = {
    local labels = ["debug", "profile"],
    local finder = {
        "OtherType": {
            "debug": "test"
        }
    },
    local make_wrapper(label, inode, pnode) = pnode + {wrapped: true},

    local my_pnode = wrapper.wrap_pnode(labels, finder, make_wrapper),
    local inode = make_inode("MyType", "test_node"),
    local result = my_pnode(inode, 1, 1),

    // Should not be wrapped because type doesn't match
    assert !std.objectHas(result, 'wrapped') : "Test 1 failed: node should not be wrapped",
    result: "PASS: No wrapping when type not in finder"
};

// Test 2: No wrapping when name substring doesn't match
local test2 = {
    local labels = ["debug"],
    local finder = {
        "MyType": {
            "debug": "nomatch"
        }
    },
    local make_wrapper(label, inode, pnode) = pnode + {wrapped: true},

    local my_pnode = wrapper.wrap_pnode(labels, finder, make_wrapper),
    local inode = make_inode("MyType", "test_node"),
    local result = my_pnode(inode, 1, 1),

    // Should not be wrapped because name doesn't contain "nomatch"
    assert !std.objectHas(result, 'wrapped') : "Test 2 failed: node should not be wrapped",
    result: "PASS: No wrapping when name doesn't match substring"
};

// Test 3: Wrapping when both type and name match
local test3 = {
    local labels = ["debug", "profile"],
    local finder = {
        "MyType": {
            "debug": "test"
        }
    },
    local make_wrapper(label, inode, pnode) = pnode + {wrapped: true, label: label},

    local my_pnode = wrapper.wrap_pnode(labels, finder, make_wrapper),
    local inode = make_inode("MyType", "test_node"),
    local result = my_pnode(inode, 1, 1),

    // Should be wrapped because type is MyType and name contains "test"
    assert std.objectHas(result, 'wrapped') : "Test 3 failed: node should be wrapped",
    assert result.label == "debug" : "Test 3 failed: wrong label",
    result: "PASS: Wrapping when type and name match"
};

// Test 4: First matching label wins
local test4 = {
    local labels = ["profile", "debug"],  // profile comes first
    local finder = {
        "MyType": {
            "debug": "test",
            "profile": "test"  // both match
        }
    },
    local make_wrapper(label, inode, pnode) = pnode + {wrapped: true, label: label},

    local my_pnode = wrapper.wrap_pnode(labels, finder, make_wrapper),
    local inode = make_inode("MyType", "test_node"),
    local result = my_pnode(inode, 1, 1),

    // Should use "profile" label because it comes first in labels array
    assert result.label == "profile" : "Test 4 failed: should use first matching label",
    result: "PASS: First matching label is used"
};

// Test 5: Label not in finder is skipped
local test5 = {
    local labels = ["notfound", "debug"],
    local finder = {
        "MyType": {
            "debug": "test"
        }
    },
    local make_wrapper(label, inode, pnode) = pnode + {wrapped: true, label: label},

    local my_pnode = wrapper.wrap_pnode(labels, finder, make_wrapper),
    local inode = make_inode("MyType", "test_node"),
    local result = my_pnode(inode, 1, 1),

    // Should use "debug" label, skipping "notfound"
    assert result.label == "debug" : "Test 5 failed: should skip missing label and use next",
    result: "PASS: Missing labels are skipped"
};

// Test 6: Wrapped pnode can be used with replace_pnode
local test6 = {
    local labels = ["debug"],
    local finder = {
        "MyType": {
            "debug": "test"
        }
    },
    local make_wrapper(label, inode, pnode) = pnode + {wrapped: true},

    local new_pg = wrapper.replace_pnode(labels, finder, make_wrapper),

    local inode = make_inode("MyType", "test_node"),
    local result = new_pg.pnode(inode, 1, 1),

    // Should be wrapped when using replaced pnode in new_pg
    assert std.objectHas(result, 'wrapped') : "Test 6 failed: replace_pnode should work",
    result: "PASS: replace_pnode integrates wrapped function"
};

// Test 7: Substring matching is case-sensitive
local test7 = {
    local labels = ["debug"],
    local finder = {
        "MyType": {
            "debug": "TEST"  // uppercase
        }
    },
    local make_wrapper(label, inode, pnode) = pnode + {wrapped: true},

    local my_pnode = wrapper.wrap_pnode(labels, finder, make_wrapper),
    local inode = make_inode("MyType", "test_node"),  // lowercase
    local result = my_pnode(inode, 1, 1),

    // Should not match because std.findSubstr is case-sensitive
    assert !std.objectHas(result, 'wrapped') : "Test 7 failed: substring match should be case-sensitive",
    result: "PASS: Substring matching is case-sensitive"
};

// Test 8: Multiple types in finder
local test8 = {
    local labels = ["debug"],
    local finder = {
        "TypeA": {
            "debug": "nodeA"
        },
        "TypeB": {
            "debug": "nodeB"
        }
    },
    local make_wrapper(label, inode, pnode) = pnode + {wrapped: true, matched_type: inode.type},

    local my_pnode = wrapper.wrap_pnode(labels, finder, make_wrapper),
    local inodeA = make_inode("TypeA", "nodeA_1"),
    local inodeB = make_inode("TypeB", "nodeB_1"),
    local resultA = my_pnode(inodeA, 1, 1),
    local resultB = my_pnode(inodeB, 1, 1),

    assert resultA.matched_type == "TypeA" : "Test 8 failed: TypeA should match",
    assert resultB.matched_type == "TypeB" : "Test 8 failed: TypeB should match",
    result: "PASS: Multiple types can be matched independently"
};

// Collect all test results
{
    tests: [
        test1.result,
        test2.result,
        test3.result,
        test4.result,
        test5.result,
        test6.result,
        test7.result,
        test8.result,
    ],
    summary: "All %d tests passed!" % std.length(self.tests)
}
