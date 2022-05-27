import React from 'react';

const AuthContext = React.createContext();
AuthContext.displayName = 'AuthContext';

function AuthProvider(props) {
  const [user, setUser] = React.useState(localStorage.user ? JSON.parse(localStorage.user) : null);

  const value = [user, setUser];
  return (
    // eslint-disable-next-line react/jsx-props-no-spreading
    <AuthContext.Provider value={value} {...props} />
  );
}

function useAuth() {
  const context = React.useContext(AuthContext);
  if (context === undefined) {
    throw new Error('useAuth must be used within a AuthProvider');
  }
  return context;
}

export { AuthProvider, useAuth };
