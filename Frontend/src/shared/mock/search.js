const DUMMY_TEAMS = [
  {
    id: "1",
    name: "Team 1 - 7K7K",
    visibility: "public",
    updated: "17th March", // Target to change to some date format
    institution: "Technische Universität München",
    membersCount: 18,
    members: [
      // not all are needed but around 3
      {
        name: "Julian",
        image: "/profileImage1.jpg",
      },
      {
        name: "Ronald",
        image: "/profileImage2.jpg",
      },
      {
        name: "Killian",
        image: "/profileImage3.jpg",
      },
    ],
  },
  {
    id: "2",
    name: "Team 2 - FF77",
    visibility: "private",
    updated: "24th March", // Target to change to some date format
    institution: "Technische Universität München",
    membersCount: 9,
    members: [
      // not all are needed but around 3
      {
        name: "Julian",
        image: "/profileImage1.jpg",
      },
      {
        name: "Ronald",
        image: "/profileImage2.jpg",
      },
      {
        name: "Killian",
        image: "/profileImage3.jpg",
      },
    ],
  },
  {
    id: "3",
    name: "Team 3 - 5623",
    visibility: "private",
    updated: "24th March", // Target to change to some date format
    institution: "Helmholtz Institut",
    membersCount: 9,
    members: [
      // not all are needed but around 3
      {
        name: "Julian",
        image: "/profileImage1.jpg",
      },
      {
        name: "Ronald",
        image: "/profileImage2.jpg",
      },
      {
        name: "Killian",
        image: "/profileImage3.jpg",
      },
    ],
  },
];

const DUMMY_INSTITUTIONS = [
  {
    id: "1",
    name: "Technische Universität München",
    logo: "/profileImage.jpg",
    updated: "17th March", // Target to change to some date format
    membersCount: 18,
    teamsCount: 20,
    members: [
      // not all are needed but around 3
      {
        name: "Julian",
        image: "/profileImage1.jpg",
      },
      {
        name: "Ronald",
        image: "/profileImage2.jpg",
      },
      {
        name: "Killian",
        image: "/profileImage3.jpg",
      },
    ],
  },
  {
    id: "2",
    name: "Helmholtz Institut",
    logo: "/profileImage.jpg",
    updated: "24th March", // Target to change to some date format
    membersCount: 18,
    teamsCount: 20,
    members: [
      // not all are needed but around 3
      {
        name: "Julian",
        image: "/profileImage1.jpg",
      },
      {
        name: "Ronald",
        image: "/profileImage2.jpg",
      },
      {
        name: "Killian",
        image: "/profileImage3.jpg",
      },
    ],
  },
  {
    id: "3",
    name: "Ludwig-Maximilians-Universität München",
    logo: "/profileImage.jpg",
    updated: "24th April", // Target to change to some date format
    membersCount: 10,
    teamsCount: 22,
    members: [
      // not all are needed but around 3
      {
        name: "Julian",
        image: "/profileImage1.jpg",
      },
      {
        name: "Ronald",
        image: "/profileImage2.jpg",
      },
      {
        name: "Killian",
        image: "/profileImage3.jpg",
      },
    ],
  },
];

async function querySearch(type) {
  switch (type) {
    case "Teams":
      return DUMMY_TEAMS;
    case "Institutions":
      return DUMMY_INSTITUTIONS;
  }
}

export default querySearch;
