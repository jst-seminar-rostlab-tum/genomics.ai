async function querySearch() {
    return [
        {
            id: "1",
            name: "Team 1 - 7K7K",
            visibility: "public",
            updated: "17th March", // Target to change to some date format
            institution: "Technische Universität München",
            membersCount: 18,
            members: [ // not all are needed but around 3
                {
                    name: "Julian"
                },
                {
                    name: "Ronald"
                },
                {
                    name: "Killian"
                }
            ]
        },
        {
            id: "2",
            name: "Team 2 - FF77",
            visibility: "private",
            updated: "24th March", // Target to change to some date format
            institution: "Technische Universität München",
            membersCount: 9,
            members: [
                {
                    name: "Julian"
                },
                {
                    name: "Ronald"
                },
                {
                    name: "Killian"
                }
            ]
        },
        {
            id: "3",
            name: "Team 3 - 5623",
            visibility: "private",
            updated: "24th March", // Target to change to some date format
            institution: "Helmholtz Institut",
            membersCount: 9,
            members: [
                {
                    name: "Julian"
                },
                {
                    name: "Ronald"
                },
                {
                    name: "Killian"
                }
            ]
        }
    ]
}

export default querySearch;