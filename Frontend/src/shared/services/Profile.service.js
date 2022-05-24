import MockProfileService from './mock/Profile.service';
import axiosInstance from './axiosInstance';

const MOCK_PROFILE = false;

let _cachedProfile;

const ProfileService = MOCK_PROFILE ? MockProfileService : {
  getProfile: async (options) => {
    const allowCache = (options || {}).allowCache ?? true;
    if (allowCache && _cachedProfile) return _cachedProfile;
    const { data } = await axiosInstance.get('/profile');
    _cachedProfile = data;
    _cachedProfile.id = _cachedProfile._id;
    return _cachedProfile;
  },

  clearProfileCache: () => {
    _cachedProfile = null;
  },
};

export default ProfileService;
