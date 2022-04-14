import { BACKEND_ADDRESS } from '../utils/common/constants';
import { getAuthAndJsonHeader } from '../utils/common/utils';

export default async function updateProfile(userData, newPassword, setUser) {
  const body = {
    first_name: userData.firstName,
    last_name: userData.lastName,
    // email: userData.emailAddress,
    note: userData.academicAffiliation,
  };
  if (!!newPassword && newPassword !== '') {
    body.password = newPassword;
  }
  return fetch(`${BACKEND_ADDRESS}/update_profile`, {
    method: 'POST',
    headers: getAuthAndJsonHeader(),
    body: JSON.stringify(body),
  }).then((response) => {
    if (response.status !== 200) {
      throw Error("Couldn't save changes");
    }
    setUser((prevState) => ({
      ...prevState,
      email: userData.emailAddress,
      firstName: userData.firstName,
      lastName: userData.lastName,
      note: userData.academicAffiliation,
    }));
  });
}
