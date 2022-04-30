import axiosInstance from './axiosInstance';

const MODEL = 'users';

const UserService = {
  getOwnTeams: async () => {
    const { data } = await axiosInstance.get(`/${MODEL}/ownteams`);
    return data;
  },
};

export default UserService;
