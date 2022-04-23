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
    logo: "	https://d3r623tes721q0.cloudfront.net/gIAPUlMwpW3uW_J-bf8e03e8590b237d198916dc62b1a43dce4ce066/favicon_profile.png",
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
    logo: "https://www.hzdr.de/db/Pic?pOid=55058",
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
    logo: "/profileImage.png",
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

const DUMMY_USERS = [
  {
    id: "1",
    name: "Max Musterman",
    email: "max.mustermann@tum.de",
    image: "/profileImage1.jpg",
    affiliation: "ELTE",
  },
  {
    id: "2",
    name: "John Smith",
    email: "john.smith@tum.de",
    image: "/profileImage2.jpg",
    affiliation: "ELTE",
  },
  {
    id: "3",
    name: "Alison Henderson",
    email: "alison.henderson@tum.de",
    image: "/profileImage3.jpg",
    affiliation: "ELTE",
  },
];

const DUMMY_PROJECTS = [
  {
    id: "1",
    name: "Job f567",
    type: "Gene Mapper",
    institution: "Helmholtz Institut",
    team: "Heinzig Lab",
  },
  {
    id: "2",
    name: "Job f568",
    type: "Gene Mapper",
    institution: "Technische Universität München",
    team: "Rost Lab",
  },
  {
    id: "3",
    name: "Job f569",
    type: "Gene Mapper",
    institution: "Ludwig-Maximilians-Universität München",
    team: "Med Lab",
  },
  {
    id: "4",
    name: "Job f709",
    type: "Gene Mapper",
    team: "Med Lab",
  },
];

async function querySearch(type, keyword) {
  let requestedData = [];
  switch (type) {
    case "teams":
      requestedData = DUMMY_TEAMS;
      break;
    case "institutions":
      requestedData = DUMMY_INSTITUTIONS;
      break;
    case "users":
      requestedData = DUMMY_USERS;
      break;
    case "projects":
      requestedData = DUMMY_PROJECTS;
      break;
  }
  return requestedData.filter((item) =>
    item.name.toLowerCase().includes(keyword)
  );
}

export default querySearch;
