const ProfileService = {
  async getProfile() {
    console.warn('You are using the mock version of getProfile(). You might want the normal one (one directory level up).');
    return {
      id: 1,
      firstName: 'Paul',
      lastName: 'Schwind',
      email: 'paul@example.com',
      avatarUrl: 'https://pschwind.de/img/favicon.png',
    };
  },
};

export default ProfileService;
