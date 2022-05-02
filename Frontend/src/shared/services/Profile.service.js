import axiosInstance from './axiosInstance';

let _cachedProfile;

const ProfileService = {
  getProfile: async (options) => {
    const allowCache = (options || {}).allowCache ?? true;
    if (allowCache && _cachedProfile) return _cachedProfile;
    const { data } = await axiosInstance.get('/profile');
    _cachedProfile = data;
    _cachedProfile.id = _cachedProfile._id;
    return _cachedProfile;
  },
};

export default ProfileService;
