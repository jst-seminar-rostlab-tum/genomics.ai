import { BACKEND_ADDRESS } from 'shared/utils/common/constants';
import { getAuthAndJsonHeader } from 'shared/utils/common/utils';

let _cachedProfile;

export default async function getProfile(options) {
  const allowCache = (options || {}).allowCache ?? true;
  if (allowCache && _cachedProfile) return _cachedProfile;
  _cachedProfile = await fetch(`${BACKEND_ADDRESS}/profile`, {
    headers: getAuthAndJsonHeader(),
  }).then((resp) => resp.json());
  _cachedProfile.id = _cachedProfile._id;
  return _cachedProfile;
}
