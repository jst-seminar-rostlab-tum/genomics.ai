import MockProfileService from './mock/Profile.service';
import axiosInstance from './axiosInstance';

const MOCK_PROFILE = false;

let _cachedProfile;

const ProfileService = MOCK_PROFILE ? MockProfileService : {
  async getProfile(options) {
    const allowCache = (options || {}).allowCache ?? true;
    if (allowCache && _cachedProfile) return _cachedProfile;
    const { data } = await axiosInstance.get('/profile');
    _cachedProfile = data;
    _cachedProfile.id = _cachedProfile._id;
    return _cachedProfile;
  },

  clearProfileCache() {
    _cachedProfile = null;
  },

  updateProfile(userData, newPassword) {
    const body = {
      first_name: userData.firstName,
      last_name: userData.lastName,
      // email: userData.emailAddress,
      note: userData.academicAffiliation,
    };
    if (!!newPassword && newPassword !== '') {
      body.password = newPassword;
    }
    return axiosInstance.post('/update_profile', body)
      .then((response) => {
        if (response.status !== 200) {
          throw Error("Couldn't save changes");
        }
        this.clearProfileCache();
        localStorage.setItem('user', JSON.stringify({
          firstName: userData.firstName,
          lastName: userData.lastName,
          email: userData.emailAddress,
          note: userData.academicAffiliation,
        }));
      });
  },
};

export default ProfileService;
