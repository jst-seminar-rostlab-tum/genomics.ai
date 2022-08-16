/**
 *  Stores the raw data needed to create a user.
 */
export interface AddUserDTO {
  firstName: string;
  lastName: string;
  email: string;
  password: string;
  note: string;
  avatarUrl?: string;
}

/**
 *  Stores the raw data to update a user.
 */
export interface UpdateUserDTO {
  firstName?: string;
  lastName?: string;
  email?: string; // TODO Changing email-address is not implemented yet.
  password?: string;
  note?: string;
  avatarUrl?: string;
}
