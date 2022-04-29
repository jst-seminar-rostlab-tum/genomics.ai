export default async function getMember(memberId) {
  return {
    id: memberId,
    firstName: ['', 'Paul', 'Julian', 'Vigan', 'Nikita', 'Cem', 'Ronald', 'Amin'][memberId],
    lastName: ['', 'Schwind', 'Giebisch', 'Lladrovci', 'Charushnikov', 'Sarıca', 'Skorobogat', 'Ben Saad'][memberId],
    email: ['', 'paul@example.com', 'julian@example.com', 'vigan@example.com', 'nikita@example.com', 'cem@example.com', 'ronald@example.com', 'amin@example.com'][memberId],
    profilePictureURL: [null, null, null, null, null, null, 'https://www.genecruncher.com/memberphotos/Frontend-RonaldSkorobogat.jpg', 'https://www.genecruncher.com/memberphotos/Frontend_Amin%20Ben%20Saad.png'][memberId],
  };
}
